// System bit-exact check: sg_top (N=2) vs system golden model.  vvp tb_top.vvp +PFX=results/rtl_vectors/sys_inj +N=86400
// Compile with -Ptb_top.K=1 or 2 (default 2).
`timescale 1ns/1ps
module tb_top #(parameter K = 2);
    integer N, i, ea0, ea1, em_p, em_b, e_err, kk;
    reg [1023:0] pfx, fn;
    reg [31:0] x0 [0:99999]; reg [31:0] x1 [0:99999]; reg v0 [0:99999]; reg v1 [0:99999];
    reg a0 [0:99999]; reg a1 [0:99999]; reg [1:0] mm [0:99999];
    reg clk = 0, rst_n = 0, tick = 0; reg [1:0] valid = 0; reg [63:0] x = 0;
    wire [1:0] al_p, al_b, mode_p, mode_b; wire ev_p, ev_b, sf_p, sf_b, pr_p, pr_b, err_p, err_b;
    sg_top #(.N(2), .K(K), .DET_TMR(1), .CTRL_TMR(1)) prot (.clk(clk), .rst_n(rst_n), .tick(tick), .valid(valid), .x(x),
        .alarm(al_p), .event_o(ev_p), .mode(mode_p), .scrub_fast(sf_p), .protect(pr_p), .err(err_p));
    sg_top #(.N(2), .K(K), .DET_TMR(0), .CTRL_TMR(0)) bare (.clk(clk), .rst_n(rst_n), .tick(tick), .valid(valid), .x(x),
        .alarm(al_b), .event_o(ev_b), .mode(mode_b), .scrub_fast(sf_b), .protect(pr_b), .err(err_b));
    always #5 clk = ~clk;
    initial begin
        if (!$value$plusargs("N=%d", N)) N = 86400;
        if (!$value$plusargs("PFX=%s", pfx)) pfx = "vec";
        $sformat(fn, "%0s_c0_x.hex", pfx); $readmemh(fn, x0);   $sformat(fn, "%0s_c1_x.hex", pfx); $readmemh(fn, x1);
        $sformat(fn, "%0s_c0_v.hex", pfx); $readmemh(fn, v0);   $sformat(fn, "%0s_c1_v.hex", pfx); $readmemh(fn, v1);
        $sformat(fn, "%0s_c0_alarm.hex", pfx); $readmemh(fn, a0); $sformat(fn, "%0s_c1_alarm.hex", pfx); $readmemh(fn, a1);
        kk = K; $sformat(fn, "%0s_mode_k%0d.hex", pfx, kk); $readmemh(fn, mm);
        ea0 = 0; ea1 = 0; em_p = 0; em_b = 0; e_err = 0;
        repeat (3) @(negedge clk); rst_n = 1;
        for (i = 0; i < N; i = i + 1) begin
            @(negedge clk); x = {x1[i], x0[i]}; valid = {v1[i], v0[i]}; tick = 1;
            @(negedge clk);
            if (al_p[0] !== a0[i] || al_b[0] !== a0[i]) ea0 = ea0 + 1;
            if (al_p[1] !== a1[i] || al_b[1] !== a1[i]) ea1 = ea1 + 1;
            if (mode_p !== mm[i]) em_p = em_p + 1;
            if (mode_b !== mm[i]) em_b = em_b + 1;
            if (err_p) e_err = e_err + 1;
            tick = 0; valid = 0;
        end
        $display("SYS %0s K=%0d N=%0d alarm_err(c0=%0d c1=%0d) mode_err(prot=%0d bare=%0d) spurious_err_flags=%0d",
                 pfx, K, N, ea0, ea1, em_p, em_b, e_err);
        $finish;
    end
endmodule
