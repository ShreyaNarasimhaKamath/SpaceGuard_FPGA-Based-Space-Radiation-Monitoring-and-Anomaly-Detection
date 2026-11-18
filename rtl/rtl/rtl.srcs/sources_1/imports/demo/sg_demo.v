// SpaceGuard ZCU104 live demo. Board independent core (sg_demo_core) + thin board top (sg_board_top).
//  * replays 9,000 real GOES-16 samples (with an injected radiation event near sample 4740) from block RAM
//  * runs a PROTECTED (TMR) detector+controller and an UNPROTECTED one side by side on the same samples
//  * a push button (or the auto switch at sample AUTO_AT) flips ONE bit of channel 0's baseline register in both
//  * LEDs: 0 protected not-NORMAL, 1 protected PROTECT, 2 unprotected not-NORMAL, 3 TMR corrected an upset (stretched)
// Demo-only fault injection lives here; the verified RTL (sg_top, sg_detector, sg_ctrl, sg_next) is not modified.
`timescale 1ns/1ps

module sg_detector_tmr_fi #(parameter FREEZE = 1, parameter [65:0] MASK = 66'd1 << 26) (
    input  wire        clk, rst_n, valid, inj,
    input  wire [31:0] x,
    output wire        alarm, tmr_err
);
    (* DONT_TOUCH = "TRUE" *) reg [65:0] st0, st1, st2;
    wire [65:0] vst = (st0 & st1) | (st1 & st2) | (st0 & st2);
    wire [65:0] n0, n1, n2;
    (* DONT_TOUCH = "TRUE" *) sg_next #(.FREEZE(FREEZE)) u0 (.st(vst), .valid(valid), .x(x), .nst(n0));
    (* DONT_TOUCH = "TRUE" *) sg_next #(.FREEZE(FREEZE)) u1 (.st(vst), .valid(valid), .x(x), .nst(n1));
    (* DONT_TOUCH = "TRUE" *) sg_next #(.FREEZE(FREEZE)) u2 (.st(vst), .valid(valid), .x(x), .nst(n2));
    always @(posedge clk) begin
        st0 <= rst_n ? n0 : 66'd0;
        st1 <= rst_n ? (n1 ^ (inj ? MASK : 66'd0)) : 66'd0;   // upset lands in ONE copy only
        st2 <= rst_n ? n2 : 66'd0;
    end
    assign alarm   = vst[64];
    assign tmr_err = (st0 != vst) | (st1 != vst) | (st2 != vst);
endmodule

module sg_detector_simplex_fi #(parameter FREEZE = 1, parameter [65:0] MASK = 66'd1 << 26) (
    input  wire        clk, rst_n, valid, inj,
    input  wire [31:0] x,
    output wire        alarm
);
    reg  [65:0] st0;
    wire [65:0] n0;
    (* DONT_TOUCH = "TRUE" *) sg_next #(.FREEZE(FREEZE)) u0 (.st(st0), .valid(valid), .x(x), .nst(n0));
    always @(posedge clk) st0 <= rst_n ? (n0 ^ (inj ? MASK : 66'd0)) : 66'd0;   // the only copy is hit
    assign alarm = st0[64];
endmodule

module sg_demo_core #(
    parameter NS        = 9000,          // samples in the ROM
    parameter AUTO_AT   = 4600,          // sample at which the auto switch injects the upset
    parameter TICK_SLOW = 125_000_000,   // clocks per sample, dip[0]=1: real time, 1 sample per second at 125 MHz
    parameter TICK_FAST = 625_000,       // dip[0]=0: 200 samples per second (a 9,000-sample loop takes 45 s)
    parameter DEB       = 2_500_000,     // button debounce, clocks (20 ms at 125 MHz)
    parameter STRETCH   = 62_500_000,    // TMR-flag LED hold, clocks (0.5 s)
    parameter ROMFILE   = "sg_demo_rom.mem"
) (
    input  wire       clk,
    input  wire [1:0] dip,               // dip[0]: 1 = real time, 0 = fast.  dip[1]: 1 = auto-inject at AUTO_AT
    input  wire       btn_inj, btn_rst,
    output wire [3:0] led,
    // debug / simulation taps
    output wire [1:0] mode_p, mode_b,
    output wire       tick_o, inj_o, tmr_flag,
    output wire [15:0] addr_o
);
    // ---- input conditioning ----
    reg [1:0] dsy0, dsy1; reg [1:0] bi, br;
    always @(posedge clk) begin dsy0 <= dip; dsy1 <= dsy0; bi <= {bi[0], btn_inj}; br <= {br[0], btn_rst}; end
    reg [31:0] dbi = 0, dbr = 0; reg inj_lvl = 0, rst_lvl = 0, inj_lvl_q = 0;
    always @(posedge clk) begin
        if (bi[1] == inj_lvl) dbi <= 0; else if (dbi == DEB) begin inj_lvl <= bi[1]; dbi <= 0; end else dbi <= dbi + 1;
        if (br[1] == rst_lvl) dbr <= 0; else if (dbr == DEB) begin rst_lvl <= br[1]; dbr <= 0; end else dbr <= dbr + 1;
        inj_lvl_q <= inj_lvl;
    end
    wire inj_btn = inj_lvl & ~inj_lvl_q;

    // ---- reset: power-on, button, and loop wrap ----
    reg [7:0] por = 0; reg [3:0] wrap_cnt = 0;
    always @(posedge clk) begin
        if (!(&por)) por <= por + 1;
    end
    wire rst_n = (&por) & ~rst_lvl & (wrap_cnt == 0);

    // ---- ROM replay ----
    (* rom_style = "block" *) reg [65:0] rom [0:NS-1];
    initial $readmemh(ROMFILE, rom);
    reg [65:0] rdata; reg [15:0] addr = 0; reg [31:0] cnt = 0;
    reg tick = 0; reg [1:0] valid = 0; reg [63:0] x = 0; reg inj_auto = 0;
    wire [31:0] period = dsy1[0] ? TICK_SLOW : TICK_FAST;
    always @(posedge clk) begin
        rdata <= rom[addr];
        tick <= 1'b0; valid <= 2'b00; inj_auto <= 1'b0;
        if (wrap_cnt != 0) wrap_cnt <= wrap_cnt - 1;
        if (!rst_n) begin cnt <= 0; addr <= 0; end
        else begin
            cnt <= (cnt >= period - 1) ? 32'd0 : cnt + 1;
            if (cnt == 2) begin
                tick <= 1'b1; valid <= {rdata[65], rdata[64]}; x <= rdata[63:0];
                if (dsy1[1] && addr == AUTO_AT) inj_auto <= 1'b1;
            end
            if (cnt == 3) begin
                if (addr == NS - 1) begin addr <= 0; wrap_cnt <= 4'd8; end else addr <= addr + 1;
            end
        end
    end
    wire inj = inj_btn | inj_auto;

    // ---- protected pipeline (TMR detectors + TMR controller) ----
    wire [1:0] ap; wire [1:0] ep; wire cerr_p;
    sg_detector_tmr_fi dp0 (.clk(clk), .rst_n(rst_n), .valid(valid[0]), .inj(inj), .x(x[31:0]),  .alarm(ap[0]), .tmr_err(ep[0]));
    sg_detector_tmr_fi dp1 (.clk(clk), .rst_n(rst_n), .valid(valid[1]), .inj(1'b0), .x(x[63:32]), .alarm(ap[1]), .tmr_err(ep[1]));
    sg_ctrl #(.TMR(1)) cp (.clk(clk), .rst_n(rst_n), .tick(tick), .ev_in(&ap), .valid_any(|valid), .mode(mode_p), .err(cerr_p));

    // ---- unprotected pipeline (single copies) ----
    wire [1:0] ab; wire cerr_b;
    sg_detector_simplex_fi db0 (.clk(clk), .rst_n(rst_n), .valid(valid[0]), .inj(inj), .x(x[31:0]),  .alarm(ab[0]));
    sg_detector_simplex_fi db1 (.clk(clk), .rst_n(rst_n), .valid(valid[1]), .inj(1'b0), .x(x[63:32]), .alarm(ab[1]));
    sg_ctrl #(.TMR(0)) cb (.clk(clk), .rst_n(rst_n), .tick(tick), .ev_in(&ab), .valid_any(|valid), .mode(mode_b), .err(cerr_b));

    // ---- LEDs ----
    wire err_any = |ep | cerr_p;
    reg [31:0] str = 0;
    always @(posedge clk) if (err_any) str <= STRETCH; else if (str != 0) str <= str - 1;
    assign led = {str != 0, mode_b != 2'd0, mode_p[1], mode_p != 2'd0};
    assign tick_o = tick; assign inj_o = inj; assign tmr_flag = err_any; assign addr_o = addr;
endmodule

// Board top for the ZCU104: 125 MHz differential system clock.
module sg_board_top (
    input  wire       clk_p, clk_n,
    input  wire [1:0] dip,
    input  wire       btn_inj, btn_rst,
    output wire [3:0] led
);
    wire clk;
    IBUFDS ib (.I(clk_p), .IB(clk_n), .O(clk));
    sg_demo_core core (.clk(clk), .dip(dip), .btn_inj(btn_inj), .btn_rst(btn_rst), .led(led),
        .mode_p(), .mode_b(), .tick_o(), .inj_o(), .tmr_flag(), .addr_o());
endmodule
