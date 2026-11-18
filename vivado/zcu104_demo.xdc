# SpaceGuard ZCU104 demo constraints.  Top module: sg_board_top
# 300 MHz system clock CLK_300_P/N (pins and standard from Vivado's ZCU104 board file part0_pins.xml v1.1).
# The design divides it by 4 (BUFGCE_DIV) to 75 MHz; Vivado propagates the clock through the divider.
set_property -dict {PACKAGE_PIN AH18 IOSTANDARD DIFF_SSTL12} [get_ports clk_p]
set_property -dict {PACKAGE_PIN AH17 IOSTANDARD DIFF_SSTL12} [get_ports clk_n]
create_clock -name clk300 -period 3.333 [get_ports clk_p]

# ---- USER I/O: pins from Vivado's ZCU104 board file (part0_pins.xml v1.1), all LVCMOS33 ----
# LEDs: led[0] protected not-NORMAL, led[1] protected PROTECT, led[2] unprotected not-NORMAL, led[3] TMR corrected an upset
set_property -dict {PACKAGE_PIN D5 IOSTANDARD LVCMOS33} [get_ports {led[0]}]   ;# GPIO_LED_0_LS
set_property -dict {PACKAGE_PIN D6 IOSTANDARD LVCMOS33} [get_ports {led[1]}]   ;# GPIO_LED_1_LS
set_property -dict {PACKAGE_PIN A5 IOSTANDARD LVCMOS33} [get_ports {led[2]}]   ;# GPIO_LED_2_LS
set_property -dict {PACKAGE_PIN B5 IOSTANDARD LVCMOS33} [get_ports {led[3]}]   ;# GPIO_LED_3_LS
# DIP switches: dip[0]=1 real time (0 = fast), dip[1]=1 auto-inject at sample 4600
set_property -dict {PACKAGE_PIN E4 IOSTANDARD LVCMOS33} [get_ports {dip[0]}]   ;# GPIO_DIP_SW0
set_property -dict {PACKAGE_PIN D4 IOSTANDARD LVCMOS33} [get_ports {dip[1]}]   ;# GPIO_DIP_SW1
# Push buttons: btn_inj = inject one upset, btn_rst = restart the replay
set_property -dict {PACKAGE_PIN B4 IOSTANDARD LVCMOS33} [get_ports btn_inj]    ;# GPIO_PB_SW0
set_property -dict {PACKAGE_PIN C4 IOSTANDARD LVCMOS33} [get_ports btn_rst]    ;# GPIO_PB_SW1
set_false_path -from [get_ports {dip[*] btn_inj btn_rst}]
set_false_path -to   [get_ports {led[*]}]
