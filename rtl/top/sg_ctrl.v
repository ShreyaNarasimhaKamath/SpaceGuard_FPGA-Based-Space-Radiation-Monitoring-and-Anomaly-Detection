// SpaceGuard mitigation controller: turns the channel-vote event into an autonomous protection level.
//   mode 0 NORMAL  : nominal scrub rate
//   mode 1 ALERT   : scrub_fast
//   mode 2 PROTECT : scrub_fast + protect (heavier protection downstream); mode 3 is illegal and treated as PROTECT
// Escalation: event persists ESC ticks -> PROTECT. De-escalation is staged (PROTECT -> ALERT -> NORMAL), each stage
// held HOLD ticks after the event clears. Fail-safe: if no channel delivers a valid sample for WD ticks, go to PROTECT.
// State (30 bits): [1:0] mode  [11:2] hold  [21:12] esc  [29:22] wd.   Time base: one 'tick' per sample period.

module sg_ctrl_next #(parameter HOLD = 600, ESC = 300, WD = 60) (
    input  wire [29:0] st,
    input  wire        tick, ev_in, valid_any,
    output reg  [29:0] nst
);
    wire [1:0] mode = st[1:0];
    wire [9:0] hold = st[11:2];
    wire [9:0] esc  = st[21:12];
    wire [7:0] wd   = st[29:22];
    wire [7:0] wd_n = valid_any ? 8'd0 : ((wd == 8'd255) ? wd : wd + 8'd1);
    wire       wdf  = (wd_n >= WD);
    wire       ev   = ev_in | wdf;
    wire [9:0] esc_n = (esc == 10'd1023) ? esc : esc + 10'd1;
    reg  [1:0] m; reg [9:0] h, e;
    always @* begin
        nst = st; m = mode; h = hold; e = esc;
        if (tick) begin
            case (mode)
                2'd0: if (ev) begin m = 2'd1; h = HOLD; e = 10'd0; end
                2'd1: if (ev) begin h = HOLD; e = esc_n; if (esc_n >= ESC) m = 2'd2; end
                      else if (hold <= 10'd1) begin m = 2'd0; h = 10'd0; e = 10'd0; end
                      else h = hold - 10'd1;
                default: if (ev) h = HOLD;
                      else if (hold <= 10'd1) begin m = 2'd1; h = HOLD; e = 10'd0; end
                      else h = hold - 10'd1;
            endcase
            if (wdf) begin m = 2'd2; h = HOLD; end
            nst = {wd_n, e, h, m};
        end
    end
endmodule

// TMR=1: three state copies, whole-state bitwise vote, voted state fed back (self-correcting). TMR=0: one copy only
// (copies 1 and 2 are not instantiated, so the unprotected variant is a fair area/timing baseline).
module sg_ctrl #(parameter TMR = 1, HOLD = 600, ESC = 300, WD = 60) (
    input  wire clk, rst_n, tick, ev_in, valid_any,
    output wire [1:0] mode,
    output wire err
);
    (* DONT_TOUCH = "TRUE" *) reg [29:0] st0;
    wire [29:0] vst, n0;
    wire        err_w;
    generate if (TMR) begin : g_t
        (* DONT_TOUCH = "TRUE" *) reg [29:0] st1, st2;
        wire [29:0] n1, n2;
        assign vst = (st0 & st1) | (st1 & st2) | (st0 & st2);
        (* DONT_TOUCH = "TRUE" *) sg_ctrl_next #(HOLD, ESC, WD) u1 (.st(vst), .tick(tick), .ev_in(ev_in), .valid_any(valid_any), .nst(n1));
        (* DONT_TOUCH = "TRUE" *) sg_ctrl_next #(HOLD, ESC, WD) u2 (.st(vst), .tick(tick), .ev_in(ev_in), .valid_any(valid_any), .nst(n2));
        always @(posedge clk) begin
            st1 <= rst_n ? n1 : 30'd0;
            st2 <= rst_n ? n2 : 30'd0;
        end
        assign err_w = (st0 != vst) | (st1 != vst) | (st2 != vst);
    end else begin : g_s
        assign vst   = st0;
        assign err_w = 1'b0;
    end endgenerate
    (* DONT_TOUCH = "TRUE" *) sg_ctrl_next #(HOLD, ESC, WD) u0 (.st(vst), .tick(tick), .ev_in(ev_in), .valid_any(valid_any), .nst(n0));
    always @(posedge clk) st0 <= rst_n ? n0 : 30'd0;
    assign mode = vst[1:0];
    assign err  = err_w;
endmodule
