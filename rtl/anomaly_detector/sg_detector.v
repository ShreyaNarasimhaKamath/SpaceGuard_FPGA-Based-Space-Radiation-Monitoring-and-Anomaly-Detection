// SpaceGuard single-channel detector, two variants sharing sg_next:
//   sg_detector_simplex : one state copy
//   sg_detector_tmr     : three state copies, bitwise majority vote on the WHOLE state every cycle,
//                         voted state feeds back to all three copies (self-correcting, errors cannot accumulate).
// Outputs: alarm (voted 'on'), base_o (voted baseline, debug), tmr_err (pulse: a copy disagreed with the vote).

module sg_detector_simplex #(parameter FREEZE = 1) (
    input  wire        clk, rst_n, valid,
    input  wire [31:0] x,
    output wire        alarm,
    output wire [31:0] base_o
);
    reg  [65:0] st0;
    wire [65:0] n0;
    (* DONT_TOUCH = "TRUE" *) sg_next #(.FREEZE(FREEZE)) u0 (.st(st0), .valid(valid), .x(x), .nst(n0));
    always @(posedge clk) st0 <= rst_n ? n0 : 66'd0;
    assign alarm  = st0[64];
    assign base_o = st0[31:0];
endmodule

module sg_detector_tmr #(parameter FREEZE = 1) (
    input  wire        clk, rst_n, valid,
    input  wire [31:0] x,
    output wire        alarm,
    output wire [31:0] base_o,
    output wire        tmr_err
);
    (* DONT_TOUCH = "TRUE" *) reg  [65:0] st0, st1, st2;
    wire [65:0] vst = (st0 & st1) | (st1 & st2) | (st0 & st2);
    wire [65:0] n0, n1, n2;
    (* DONT_TOUCH = "TRUE" *) sg_next #(.FREEZE(FREEZE)) u0 (.st(vst), .valid(valid), .x(x), .nst(n0));
    (* DONT_TOUCH = "TRUE" *) sg_next #(.FREEZE(FREEZE)) u1 (.st(vst), .valid(valid), .x(x), .nst(n1));
    (* DONT_TOUCH = "TRUE" *) sg_next #(.FREEZE(FREEZE)) u2 (.st(vst), .valid(valid), .x(x), .nst(n2));
    always @(posedge clk) begin
        st0 <= rst_n ? n0 : 66'd0;
        st1 <= rst_n ? n1 : 66'd0;
        st2 <= rst_n ? n2 : 66'd0;
    end
    assign alarm   = vst[64];
    assign base_o  = vst[31:0];
    assign tmr_err = (st0 != vst) | (st1 != vst) | (st2 != vst);
endmodule
