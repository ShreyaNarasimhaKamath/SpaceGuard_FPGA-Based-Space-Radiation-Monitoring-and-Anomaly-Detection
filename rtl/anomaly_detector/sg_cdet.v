// SpaceGuard count-channel detector (hazard channels, e.g. >2 MeV electron integral). Board independent Verilog-2001.
// Input c: integer counts per sample period. Leaky accumulator with FR=8 fractional bits (shift-only):
//   L <- L - (L>>SH) + (c<<8)   (saturating 32 bit);  alarm Schmitt: set when L >= TH, clear when L < TL.
// Defaults: SH=8 (about 256 s window), TH/TL for the 1000 pfu SWPC >2 MeV level (Q=1326 pfu per count/s): TH=49424, TL=34597.
// State (33 bits): [31:0] L, [32] on.
module sg_cdet_next #(parameter SH = 8, TH = 49424, TL = 34597) (
    input  wire [32:0] st,
    input  wire        valid,
    input  wire [15:0] c,
    output reg  [32:0] nst
);
    wire [31:0] L = st[31:0];
    wire        on = st[32];
    wire [33:0] ln = {2'b00, L} - {2'b00, (L >> SH)} + {10'd0, c, 8'd0};
    wire [31:0] L_n = (ln > 34'h0FFFFFFFF) ? 32'hFFFFFFFF : ln[31:0];
    wire        on_n = on ? (L_n >= TL) : (L_n >= TH);
    always @* begin
        nst = st;
        if (valid) nst = {on_n, L_n};
    end
endmodule

module sg_cdet_simplex #(parameter SH = 8, TH = 49424, TL = 34597) (
    input wire clk, rst_n, valid, input wire [15:0] c, output wire alarm
);
    reg [32:0] st0; wire [32:0] n0;
    (* DONT_TOUCH = "TRUE" *) sg_cdet_next #(SH, TH, TL) u0 (.st(st0), .valid(valid), .c(c), .nst(n0));
    always @(posedge clk) st0 <= rst_n ? n0 : 33'd0;
    assign alarm = st0[32];
endmodule

module sg_cdet_tmr #(parameter SH = 8, TH = 49424, TL = 34597) (
    input wire clk, rst_n, valid, input wire [15:0] c, output wire alarm, output wire tmr_err
);
    (* DONT_TOUCH = "TRUE" *) reg [32:0] st0, st1, st2;
    wire [32:0] vst = (st0 & st1) | (st1 & st2) | (st0 & st2);
    wire [32:0] n0, n1, n2;
    (* DONT_TOUCH = "TRUE" *) sg_cdet_next #(SH, TH, TL) u0 (.st(vst), .valid(valid), .c(c), .nst(n0));
    (* DONT_TOUCH = "TRUE" *) sg_cdet_next #(SH, TH, TL) u1 (.st(vst), .valid(valid), .c(c), .nst(n1));
    (* DONT_TOUCH = "TRUE" *) sg_cdet_next #(SH, TH, TL) u2 (.st(vst), .valid(valid), .c(c), .nst(n2));
    always @(posedge clk) begin
        st0 <= rst_n ? n0 : 33'd0; st1 <= rst_n ? n1 : 33'd0; st2 <= rst_n ? n2 : 33'd0;
    end
    assign alarm = vst[32];
    assign tmr_err = (st0 != vst) | (st1 != vst) | (st2 != vst);
endmodule
