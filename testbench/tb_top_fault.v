// System fault injection: one SEU every PERIOD samples into a random detector channel or the controller,
// fully protected top (TMR detectors + TMR controller) vs bare top, scored against the fault-free golden mode.
// vvp tb_top_fault.vvp +PFX=... +PERIOD=500 +SEED=1 +MODE=0   (MODE 1: same bit in two copies of the hit block, protected top)
`timescale 1ns/1ps
module tb_top_fault #(parameter K = 2);
    integer N, PERIOD, MODE, SEED, i, kk, ev, tgt, nb, nc, pb, bb, p_flag, b_flag, p_ev, b_ev, p_cyc, b_cyc, det;
    reg [1023:0] pfx, fn;
    reg [31:0] x0 [0:99999]; reg [31:0] x1 [0:99999]; reg v0 [0:99999]; reg v1 [0:99999]; reg [1:0] mm [0:99999];
    reg clk = 0, rst_n = 0, tick = 0; reg [1:0] valid = 0; reg [63:0] x = 0;
    reg [65:0] m66; reg [29:0] m30;
    wire [1:0] al_p, al_b, mode_p, mode_b; wire ev_p, ev_b, sf_p, sf_b, pr_p, pr_b, err_p, err_b;
    sg_top #(.N(2), .K(K), .DET_TMR(1), .CTRL_TMR(1)) prot (.clk(clk), .rst_n(rst_n), .tick(tick), .valid(valid), .x(x),
        .alarm(al_p), .event_o(ev_p), .mode(mode_p), .scrub_fast(sf_p), .protect(pr_p), .err(err_p));
    sg_top #(.N(2), .K(K), .DET_TMR(0), .CTRL_TMR(0)) bare (.clk(clk), .rst_n(rst_n), .tick(tick), .valid(valid), .x(x),
        .alarm(al_b), .event_o(ev_b), .mode(mode_b), .scrub_fast(sf_b), .protect(pr_b), .err(err_b));
    always #5 clk = ~clk;
    initial begin
        if (!$value$plusargs("N=%d", N)) N = 86400;
        if (!$value$plusargs("PERIOD=%d", PERIOD)) PERIOD = 500;
        if (!$value$plusargs("MODE=%d", MODE)) MODE = 0;
        if (!$value$plusargs("SEED=%d", SEED)) SEED = 1;
        if (!$value$plusargs("PFX=%s", pfx)) pfx = "vec";
        $sformat(fn, "%0s_c0_x.hex", pfx); $readmemh(fn, x0);   $sformat(fn, "%0s_c1_x.hex", pfx); $readmemh(fn, x1);
        $sformat(fn, "%0s_c0_v.hex", pfx); $readmemh(fn, v0);   $sformat(fn, "%0s_c1_v.hex", pfx); $readmemh(fn, v1);
        kk = K; $sformat(fn, "%0s_mode_k%0d.hex", pfx, kk); $readmemh(fn, mm);
        ev = 0; p_ev = 0; b_ev = 0; p_cyc = 0; b_cyc = 0; p_flag = 0; b_flag = 0; det = 0;
        repeat (3) @(negedge clk); rst_n = 1;
        for (i = 0; i < N; i = i + 1) begin
            @(negedge clk); x = {x1[i], x0[i]}; valid = {v1[i], v0[i]}; tick = 1;
            if (i > 1900 && (i % PERIOD) == 0) begin
                if (ev > 0) begin p_ev = p_ev + p_flag; b_ev = b_ev + b_flag; end
                p_flag = 0; b_flag = 0; ev = ev + 1;
                tgt = {$random(SEED)} % 3;                 // 0,1: detector channel, 2: controller
                nc  = {$random(SEED)} % 3;
                if (tgt == 2) begin
                    nb = {$random(SEED)} % 30; m30 = 30'd1 << nb;
                    bare.u_ctrl.st0 = bare.u_ctrl.st0 ^ m30;
                    if (MODE == 1) begin prot.u_ctrl.st0 = prot.u_ctrl.st0 ^ m30; prot.u_ctrl.g_t.st1 = prot.u_ctrl.g_t.st1 ^ m30; end
                    else case (nc) 0: prot.u_ctrl.st0 = prot.u_ctrl.st0 ^ m30; 1: prot.u_ctrl.g_t.st1 = prot.u_ctrl.g_t.st1 ^ m30;
                                   default: prot.u_ctrl.g_t.st2 = prot.u_ctrl.g_t.st2 ^ m30; endcase
                end else begin
                    nb = {$random(SEED)} % 66; m66 = 66'd1 << nb;
                    if (tgt == 0) begin
                        bare.g_ch[0].g_k0.g_s.d.st0 = bare.g_ch[0].g_k0.g_s.d.st0 ^ m66;
                        if (MODE == 1) begin prot.g_ch[0].g_k0.g_t.d.st0 = prot.g_ch[0].g_k0.g_t.d.st0 ^ m66; prot.g_ch[0].g_k0.g_t.d.st1 = prot.g_ch[0].g_k0.g_t.d.st1 ^ m66; end
                        else case (nc) 0: prot.g_ch[0].g_k0.g_t.d.st0 = prot.g_ch[0].g_k0.g_t.d.st0 ^ m66; 1: prot.g_ch[0].g_k0.g_t.d.st1 = prot.g_ch[0].g_k0.g_t.d.st1 ^ m66;
                                       default: prot.g_ch[0].g_k0.g_t.d.st2 = prot.g_ch[0].g_k0.g_t.d.st2 ^ m66; endcase
                    end else begin
                        bare.g_ch[1].g_k0.g_s.d.st0 = bare.g_ch[1].g_k0.g_s.d.st0 ^ m66;
                        if (MODE == 1) begin prot.g_ch[1].g_k0.g_t.d.st0 = prot.g_ch[1].g_k0.g_t.d.st0 ^ m66; prot.g_ch[1].g_k0.g_t.d.st1 = prot.g_ch[1].g_k0.g_t.d.st1 ^ m66; end
                        else case (nc) 0: prot.g_ch[1].g_k0.g_t.d.st0 = prot.g_ch[1].g_k0.g_t.d.st0 ^ m66; 1: prot.g_ch[1].g_k0.g_t.d.st1 = prot.g_ch[1].g_k0.g_t.d.st1 ^ m66;
                                       default: prot.g_ch[1].g_k0.g_t.d.st2 = prot.g_ch[1].g_k0.g_t.d.st2 ^ m66; endcase
                    end
                end
                #1; if (err_p) det = det + 1;
            end
            @(negedge clk);
            if (mode_p !== mm[i]) begin p_cyc = p_cyc + 1; p_flag = 1; end
            if (mode_b !== mm[i]) begin b_cyc = b_cyc + 1; b_flag = 1; end
            tick = 0; valid = 0;
        end
        p_ev = p_ev + p_flag; b_ev = b_ev + b_flag;
        $display("SYSFAULT %0s K=%0d mode=%0d events=%0d | bare: events_with_mode_error=%0d wrong_mode_cycles=%0d | protected: events_with_mode_error=%0d wrong_mode_cycles=%0d err_flags=%0d",
                 pfx, K, MODE, ev, b_ev, b_cyc, p_ev, p_cyc, det);
        $finish;
    end
endmodule
