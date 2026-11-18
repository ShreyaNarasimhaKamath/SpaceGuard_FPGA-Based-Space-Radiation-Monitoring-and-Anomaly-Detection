# Out-of-context clock for resource/timing estimates (300 MHz target). Not a board constraint.
create_clock -period 3.333 -name clk [get_ports clk]
