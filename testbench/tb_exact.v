// Bit-exact check: RTL (simplex and TMR) vs golden model vectors. Usage: vvp tb_exact.vvp +N=86400 +PFX=path/name
`timescale 1ns/1ps
module tb_exact;
    integer N, i, err_a_s, err_b_s, err_a_t, err_b_t, err_tmr;
    reg [1023:0] pfx; reg [1023:0] fn;
    reg [31:0] xm [0:99999]; reg vm [0:99999]; reg am [0:99999]; reg [31:0] bm [0:99999];
    reg clk = 0, rst_n = 0, valid = 0; reg [31:0] x = 0;
    wire a_s, a_t, e_t; wire [31:0] b_s, b_t;
    sg_detector_simplex ds (.clk(clk), .rst_n(rst_n), .valid(valid), .x(x), .alarm(a_s), .base_o(b_s));
    sg_detector_tmr     dt (.clk(clk), .rst_n(rst_n), .valid(valid), .x(x), .alarm(a_t), .base_o(b_t), .tmr_err(e_t));
    always #5 clk = ~clk;
    initial begin
        if (!$value$plusargs("N=%d", N)) N = 86400;
        if (!$value$plusargs("PFX=%s", pfx)) pfx = "vec";
        $sformat(fn, "%0s_x.hex", pfx);     $readmemh(fn, xm);
        $sformat(fn, "%0s_v.hex", pfx);     $readmemh(fn, vm);
        $sformat(fn, "%0s_alarm.hex", pfx); $readmemh(fn, am);
        $sformat(fn, "%0s_base.hex", pfx);  $readmemh(fn, bm);
        err_a_s = 0; err_b_s = 0; err_a_t = 0; err_b_t = 0; err_tmr = 0;
        repeat (3) @(negedge clk); rst_n = 1;
        for (i = 0; i < N; i = i + 1) begin
            @(negedge clk); x = xm[i]; valid = vm[i];
            @(negedge clk);                       // state updated at the posedge in between
            if (a_s !== am[i]) err_a_s = err_a_s + 1;
            if (b_s !== bm[i]) err_b_s = err_b_s + 1;
            if (a_t !== am[i]) err_a_t = err_a_t + 1;
            if (b_t !== bm[i]) err_b_t = err_b_t + 1;
            if (e_t) err_tmr = err_tmr + 1;
            valid = 0;
        end
        $display("EXACT %0s N=%0d simplex(alarm_err=%0d base_err=%0d) tmr(alarm_err=%0d base_err=%0d spurious_tmr_err=%0d)",
                 pfx, N, err_a_s, err_b_s, err_a_t, err_b_t, err_tmr);
        $finish;
    end
endmodule
