// SpaceGuard single-channel detector: combinational next-state function (board independent, plain Verilog-2001).
// State layout (66 bits, LSB first):
//   [31:0] base  [42:32] frz  [48:43] off  [52:49] run  [63:53] cnt  [64] on  [65] init
// x is an unsigned 32-bit fixed-point sample (F fractional bits chosen by the system, 8 in the reference model).
module sg_next #(
    parameter NU = 10,      // EMA shift when x > base (slow rise)
    parameter ND = 7,       // EMA shift when x <= base (fast fall)
    parameter M = 10,       // persistence: consecutive samples above threshold to raise alarm
    parameter R = 30,       // release: consecutive samples below threshold to clear alarm
    parameter WARM = 1800,  // valid samples before the detector is armed
    parameter TMAX = 1800,  // max samples the baseline is frozen per alarm
    parameter FREEZE = 1
) (
    input  wire [65:0] st,
    input  wire        valid,
    input  wire [31:0] x,
    output reg  [65:0] nst
);
    wire [31:0] base = st[31:0];
    wire [10:0] frz  = st[42:32];
    wire [5:0]  off  = st[48:43];
    wire [3:0]  run  = st[52:49];
    wire [10:0] cnt  = st[63:53];
    wire        on   = st[64];
    wire        init = st[65];

    // threshold = base * 1.5 using shift-add
    wire [33:0] thr  = {2'b00, base} + {3'b000, base[31:1]};
    wire        cond = (cnt >= WARM) && ({2'b00, x} > thr);

    wire [3:0]  run_n = cond ? ((run == M) ? run : run + 4'd1) : 4'd0;
    wire [5:0]  off_n = cond ? 6'd0 : ((off == R) ? off : off + 6'd1);
    wire        on_set = (!on) && (run_n >= M);
    wire        on_clr = on && (off_n >= R);
    wire        on_n   = on_set ? 1'b1 : (on_clr ? 1'b0 : on);
    wire [10:0] frz_e  = on_set ? 11'd0 : frz;
    wire        hold   = (FREEZE != 0) && on_n && (frz_e < TMAX);

    // asymmetric EMA step (arithmetic shift = floor)
    wire signed [33:0] d  = $signed({2'b00, x}) - $signed({2'b00, base});
    wire signed [33:0] ds = (d > 0) ? (d >>> NU) : (d >>> ND);
    wire signed [34:0] bn = $signed({3'b000, base}) + ds;
    wire [31:0] base_n = (bn < 0) ? 32'd0 : (bn > 35'sh0FFFFFFFF) ? 32'hFFFFFFFF : bn[31:0];

    always @* begin
        nst = st;
        if (valid) begin
            if (!init) begin
                nst = {1'b1, 1'b0, 11'd1, 4'd0, 6'd0, 11'd0, x};  // {init,on,cnt,run,off,frz,base}
            end else begin
                nst[52:49] = run_n;
                nst[48:43] = off_n;
                nst[64]    = on_n;
                nst[42:32] = hold ? frz_e + 11'd1 : frz_e;
                nst[31:0]  = hold ? base : base_n;
                nst[63:53] = (cnt == WARM) ? cnt : cnt + 11'd1;
            end
        end
    end
endmodule
