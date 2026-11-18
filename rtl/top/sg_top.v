// SpaceGuard top: N detector channels -> k-of-N channel vote -> mitigation controller. Board independent.
module sg_top #(parameter N = 2, K = 2, KIND = 0, DET_TMR = 1, CTRL_TMR = 1, HOLD = 600, ESC = 300, WD = 60) (
    input  wire           clk, rst_n,
    input  wire           tick,                 // one pulse per sample period
    input  wire [N-1:0]   valid,
    input  wire [32*N-1:0] x,
    output wire [N-1:0]   alarm,
    output wire           event_o,
    output wire [1:0]     mode,
    output wire           scrub_fast, protect,
    output wire           err                   // any TMR copy disagreed (upset seen and, if single, corrected)
);
    wire [N-1:0] terr;
    genvar c;
    generate for (c = 0; c < N; c = c + 1) begin : g_ch
        if (KIND == 0) begin : g_k0
          if (DET_TMR) begin : g_t
            sg_detector_tmr d (.clk(clk), .rst_n(rst_n), .valid(valid[c]), .x(x[32*c +: 32]),
                               .alarm(alarm[c]), .base_o(), .tmr_err(terr[c]));
          end else begin : g_s
            sg_detector_simplex d (.clk(clk), .rst_n(rst_n), .valid(valid[c]), .x(x[32*c +: 32]),
                                   .alarm(alarm[c]), .base_o());
            assign terr[c] = 1'b0;
          end
        end else begin : g_k1                       // KIND 1: count detector, x[32c+:16] = counts per sample
          if (DET_TMR) begin : g_ct
            sg_cdet_tmr d (.clk(clk), .rst_n(rst_n), .valid(valid[c]), .c(x[32*c +: 16]), .alarm(alarm[c]), .tmr_err(terr[c]));
          end else begin : g_cs
            sg_cdet_simplex d (.clk(clk), .rst_n(rst_n), .valid(valid[c]), .c(x[32*c +: 16]), .alarm(alarm[c]));
            assign terr[c] = 1'b0;
          end
        end
    end endgenerate

    integer j; reg [7:0] cnt;
    always @* begin cnt = 8'd0; for (j = 0; j < N; j = j + 1) cnt = cnt + alarm[j]; end
    assign event_o = (cnt >= K);

    wire cerr;
    sg_ctrl #(.TMR(CTRL_TMR), .HOLD(HOLD), .ESC(ESC), .WD(WD)) u_ctrl
        (.clk(clk), .rst_n(rst_n), .tick(tick), .ev_in(event_o), .valid_any(|valid), .mode(mode), .err(cerr));
    assign scrub_fast = (mode != 2'd0);
    assign protect    = mode[1];
    assign err        = (|terr) | cerr;
endmodule
