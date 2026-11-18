# SpaceGuard ZCU104 demo constraints.  Top module: sg_board_top
# 125 MHz system clock: pins taken from a community ZCU104 example (F23/E23, LVDS). VERIFY against the board part file.
set_property PACKAGE_PIN F23 [get_ports clk_p]
set_property PACKAGE_PIN E23 [get_ports clk_n]
set_property IOSTANDARD LVDS [get_ports {clk_p clk_n}]
create_clock -name clk125 -period 8.000 [get_ports clk_p]

# ---- USER I/O: PACKAGE_PINs deliberately NOT guessed. Run vivado/find_zcu104_pins.tcl, copy the LOC values here. ----
# LEDs (GPIO_LED_0..3)  led[0] protected not-NORMAL, led[1] protected PROTECT, led[2] unprotected not-NORMAL, led[3] TMR corrected
#set_property PACKAGE_PIN <GPIO_LED_0> [get_ports {led[0]}]
#set_property PACKAGE_PIN <GPIO_LED_1> [get_ports {led[1]}]
#set_property PACKAGE_PIN <GPIO_LED_2> [get_ports {led[2]}]
#set_property PACKAGE_PIN <GPIO_LED_3> [get_ports {led[3]}]
# DIP switches (GPIO_DIP_SW0..1)  dip[0]=1 real time, dip[1]=1 auto-inject at sample 4600
#set_property PACKAGE_PIN <GPIO_DIP_SW0> [get_ports {dip[0]}]
#set_property PACKAGE_PIN <GPIO_DIP_SW1> [get_ports {dip[1]}]
# Push buttons: btn_inj = inject one upset, btn_rst = restart the replay
#set_property PACKAGE_PIN <GPIO_SW_?> [get_ports btn_inj]
#set_property PACKAGE_PIN <GPIO_SW_?> [get_ports btn_rst]
# IOSTANDARD for all user I/O (3.3 V bank on the ZCU104 user I/O, confirm in the pin file): LVCMOS33
#set_property IOSTANDARD LVCMOS33 [get_ports {led[*] dip[*] btn_inj btn_rst}]
set_false_path -from [get_ports {dip[*] btn_inj btn_rst}]
set_false_path -to   [get_ports {led[*]}]
