// Checks sg_demo_core against the system golden model with the real ROM. Fast periods (8 clocks/sample) so it runs quickly.
//   run 1 (AUTO=0): no injection -> protected and unprotected mode must equal the golden mode on every sample
//   run 2 (AUTO=1): upset at sample 4600 -> protected must STILL equal golden, unprotected must differ, TMR flag must fire
`timescale 1ns/1ps
module tb_demo #(parameter AUTO = 0, parameter GOLD = "results/rtl_vectors/sys_inj_mode_k2.hex",
                 parameter ROM = "rtl/demo/sg_demo_rom.mem");
    reg clk = 0; reg [1:0] dip; wire [3:0] led; wire [1:0] mp, mb; wire tk, inj, tf; wire [15:0] ad;
    reg [1:0] gold [0:99999];
    sg_demo_core #(.TICK_SLOW(40), .TICK_FAST(8), .DEB(3), .STRETCH(50), .ROMFILE(ROM))
        dut (.clk(clk), .dip(dip), .btn_inj(1'b0), .btn_rst(1'b0), .led(led), .mode_p(mp), .mode_b(mb), .tick_o(tk), .inj_o(inj), .tmr_flag(tf), .addr_o(ad));
    always #4 clk = ~clk;
    integer s, ep, eb, nflag, ninj, last_s; reg seen;
    initial begin
        $readmemh(GOLD, gold);
        dip = {AUTO[0], 1'b0}; ep = 0; eb = 0; nflag = 0; ninj = 0; seen = 0; last_s = -1;
        // sample at the clock after tick: sample index = ad (addr still = index until cnt==3)
        repeat (9000 * 8 + 2000) begin
            @(posedge clk); #1;
            if (inj) ninj = ninj + 1;
            if (tf) nflag = nflag + 1;
            // 4 clocks after the tick the controller has settled: compare once per sample at addr's cnt==3 point
            if (tk) begin seen = 1; s = ad; end
            else if (seen) begin
                seen = 0;   // first clock after the tick pulse: controller state already updated for sample s
                if (s < 9000) begin
                    if (mp !== gold[s]) ep = ep + 1;
                    if (mb !== gold[s]) eb = eb + 1;
                end
            end
        end
        $display("DEMO auto=%0d  protected_mode_errors=%0d  unprotected_mode_errors=%0d  inj_pulses=%0d  tmr_flag_cycles=%0d", AUTO, ep, eb, ninj, nflag);
        $finish;
    end
endmodule
