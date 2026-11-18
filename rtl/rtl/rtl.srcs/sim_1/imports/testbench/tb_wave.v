// SpaceGuard waveform testbench (Vivado xsim or Icarus). Replays 9,000 real samples (with injected events) through
// the PROTECTED (TMR) and the UNPROTECTED top-level in parallel, and flips one bit of channel 0's baseline register
// in ONE copy at sample UPSET_AT. Look at the waveform: the protected alarm/mode stay correct, the bare ones go wrong.
// Vivado: add this file under Simulation Sources, set tb_wave as simulation top, check the PFX paths below, Run Simulation.
// Plain Verilog-2001. Needs the six .v files under rtl/ and the vector files in results/rtl_vectors/.
`timescale 1ns/1ps
module tb_wave #(
    parameter X0 = "C:/SpaceGuard/SpaceGuard_FPGA-Based-Space-Radiation-Monitoring-and-Anomaly-Detection/results/rtl_vectors/sys_inj_c0_x.hex",
    parameter X1 = "C:/SpaceGuard/SpaceGuard_FPGA-Based-Space-Radiation-Monitoring-and-Anomaly-Detection/results/rtl_vectors/sys_inj_c1_x.hex",
    parameter V0 = "C:/SpaceGuard/SpaceGuard_FPGA-Based-Space-Radiation-Monitoring-and-Anomaly-Detection/results/rtl_vectors/sys_inj_c0_v.hex",
    parameter V1 = "C:/SpaceGuard/SpaceGuard_FPGA-Based-Space-Radiation-Monitoring-and-Anomaly-Detection/results/rtl_vectors/sys_inj_c1_v.hex",
    parameter NS = 9000,            // samples to replay (the detector needs 1,800 to warm up)
    parameter UPSET_AT = 4600,      // sample index of the upset
    parameter UPSET_BIT = 26        // baseline bit to flip (weight 2^26 vs a baseline of about 6.5e6 = 2^22.6)
);
    reg [31:0] x0m [0:86399]; reg [31:0] x1m [0:86399]; reg v0m [0:86399]; reg v1m [0:86399];
    reg clk = 0, rst_n = 0, tick = 0; reg [1:0] valid = 0; reg [63:0] x = 0;
    reg [31:0] sample = 0; reg upset = 0;
    wire [1:0] al_p, al_b, mode_p, mode_b; wire ev_p, ev_b, sf_p, sf_b, pr_p, pr_b, err_p, err_b;

    sg_top #(.N(2), .K(2), .KIND(0), .DET_TMR(1), .CTRL_TMR(1)) prot (.clk(clk), .rst_n(rst_n), .tick(tick), .valid(valid), .x(x),
        .alarm(al_p), .event_o(ev_p), .mode(mode_p), .scrub_fast(sf_p), .protect(pr_p), .err(err_p));
    sg_top #(.N(2), .K(2), .KIND(0), .DET_TMR(0), .CTRL_TMR(0)) bare (.clk(clk), .rst_n(rst_n), .tick(tick), .valid(valid), .x(x),
        .alarm(al_b), .event_o(ev_b), .mode(mode_b), .scrub_fast(sf_b), .protect(pr_b), .err(err_b));

    // probes for the waveform window: channel-0 baseline (protected = voted value, bare = the single copy)
    wire [31:0] base_prot = prot.g_ch[0].g_k0.g_t.d.vst[31:0];
    wire [31:0] base_bare = bare.g_ch[0].g_k0.g_s.d.st0[31:0];
    wire [31:0] flux_ch0  = x[31:0];

    always #5 clk = ~clk;
    // the TMR error flag is a one-cycle pulse (the voter fixes the copy on the next edge): latch it so it shows per sample
    reg err_seen = 0;
    always @(posedge clk) if (err_p) err_seen <= 1'b1;

`ifdef CSV_OUT
    integer fcsv;
`endif
    integer i;
    initial begin
        $readmemh(X0, x0m); $readmemh(X1, x1m); $readmemh(V0, v0m); $readmemh(V1, v1m);
`ifdef CSV_OUT
        fcsv = $fopen("wave.csv", "w"); $fwrite(fcsv, "i,x0,base_prot,base_bare,alarm_prot,alarm_bare,mode_prot,mode_bare,upset,tmr_err\n");
`endif
        repeat (3) @(negedge clk); rst_n = 1;
        for (i = 0; i < NS; i = i + 1) begin
            @(negedge clk); sample = i; x = {x1m[i], x0m[i]}; valid = {v1m[i], v0m[i]}; tick = 1;
            if (i == UPSET_AT) begin
                upset = 1;
                prot.g_ch[0].g_k0.g_t.d.st1 = prot.g_ch[0].g_k0.g_t.d.st1 ^ (66'd1 << UPSET_BIT);   // one TMR copy hit
                bare.g_ch[0].g_k0.g_s.d.st0 = bare.g_ch[0].g_k0.g_s.d.st0 ^ (66'd1 << UPSET_BIT);   // the only copy hit
            end
            @(negedge clk); tick = 0; valid = 0;
`ifdef CSV_OUT
            $fwrite(fcsv, "%0d,%0d,%0d,%0d,%0d,%0d,%0d,%0d,%0d,%0d\n", i, x0m[i], base_prot, base_bare, al_p[0], al_b[0], mode_p, mode_b, upset, err_seen);
            err_seen = 0;
`endif
        end
`ifdef CSV_OUT
        $fclose(fcsv);
`endif
        $display("tb_wave done: %0d samples, upset at %0d bit %0d", NS, UPSET_AT, UPSET_BIT);
        $finish;
    end
endmodule
