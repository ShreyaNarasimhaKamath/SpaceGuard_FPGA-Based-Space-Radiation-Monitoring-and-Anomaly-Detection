// Fault-injection: single event upsets into the detector state, simplex vs TMR, scored against the fault-free golden alarm.
// vvp tb_fault.vvp +PFX=path/name +N=86400 +PERIOD=2000 +MODE=0 +SEED=1
//  MODE 0: one random state bit flipped in one random copy (TMR) / in the only copy (simplex)
//  MODE 1: the SAME bit flipped in two TMR copies at once (worst case: defeats the voter; simplex gets one flip)
`timescale 1ns/1ps
module tb_fault;
    integer N, PERIOD, MODE, SEED, i, k, nb, nc;
    integer ev, s_bad_ev, t_bad_ev, s_bad_cyc, t_bad_cyc, s_flag, t_flag, t_det;
    reg [1023:0] pfx; reg [1023:0] fn;
    reg [31:0] xm [0:99999]; reg vm [0:99999]; reg am [0:99999];
    reg clk = 0, rst_n = 0, valid = 0; reg [31:0] x = 0;
    wire a_s, a_t, e_t; wire [31:0] b_s, b_t;
    reg [65:0] mask;
    sg_detector_simplex ds (.clk(clk), .rst_n(rst_n), .valid(valid), .x(x), .alarm(a_s), .base_o(b_s));
    sg_detector_tmr     dt (.clk(clk), .rst_n(rst_n), .valid(valid), .x(x), .alarm(a_t), .base_o(b_t), .tmr_err(e_t));
    always #5 clk = ~clk;
    initial begin
        if (!$value$plusargs("N=%d", N)) N = 86400;
        if (!$value$plusargs("PERIOD=%d", PERIOD)) PERIOD = 2000;
        if (!$value$plusargs("MODE=%d", MODE)) MODE = 0;
        if (!$value$plusargs("SEED=%d", SEED)) SEED = 1;
        if (!$value$plusargs("PFX=%s", pfx)) pfx = "vec";
        $sformat(fn, "%0s_x.hex", pfx);     $readmemh(fn, xm);
        $sformat(fn, "%0s_v.hex", pfx);     $readmemh(fn, vm);
        $sformat(fn, "%0s_alarm.hex", pfx); $readmemh(fn, am);
        ev = 0; s_bad_ev = 0; t_bad_ev = 0; s_bad_cyc = 0; t_bad_cyc = 0; s_flag = 0; t_flag = 0; t_det = 0;
        repeat (3) @(negedge clk); rst_n = 1;
        for (i = 0; i < N; i = i + 1) begin
            @(negedge clk); x = xm[i]; valid = vm[i];
            if (i > 1900 && (i % PERIOD) == 0) begin           // after warm-up, one upset every PERIOD samples
                if (ev > 0) begin s_bad_ev = s_bad_ev + s_flag; t_bad_ev = t_bad_ev + t_flag; end
                s_flag = 0; t_flag = 0; ev = ev + 1;
                nb = {$random(SEED)} % 66; nc = {$random(SEED)} % 3;
                mask = 66'd1 << nb;
                ds.st0 = ds.st0 ^ mask;
                if (MODE == 1) begin dt.st0 = dt.st0 ^ mask; dt.st1 = dt.st1 ^ mask; end
                else case (nc) 0: dt.st0 = dt.st0 ^ mask; 1: dt.st1 = dt.st1 ^ mask; default: dt.st2 = dt.st2 ^ mask; endcase
                #1; if (e_t) t_det = t_det + 1;      // flag is combinational: visible before the voter overwrites the bad copy
            end
            @(negedge clk);
            if (a_s !== am[i]) begin s_bad_cyc = s_bad_cyc + 1; s_flag = 1; end
            if (a_t !== am[i]) begin t_bad_cyc = t_bad_cyc + 1; t_flag = 1; end
            valid = 0;
        end
        s_bad_ev = s_bad_ev + s_flag; t_bad_ev = t_bad_ev + t_flag;
        $display("FAULT mode=%0d period=%0d events=%0d | simplex: events_with_alarm_error=%0d wrong_alarm_cycles=%0d | tmr: events_with_alarm_error=%0d wrong_alarm_cycles=%0d tmr_err_pulses=%0d",
                 MODE, PERIOD, ev, s_bad_ev, s_bad_cyc, t_bad_ev, t_bad_cyc, t_det);
        $finish;
    end
endmodule
