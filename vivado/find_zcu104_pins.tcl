# Prints the ZCU104 user LED / DIP / button pin assignments from Vivado's own board files.
# Run in the Vivado Tcl console:  source C:/SpaceGuard/SpaceGuard_FPGA-Based-Space-Radiation-Monitoring-and-Anomaly-Detection/vivado/find_zcu104_pins.tcl
set files [glob -nocomplain C:/AMDDesignTools/2025.2/data/xhub/boards/XilinxBoardStore/boards/Xilinx/zcu104/*/part0_pins.xml]
if {[llength $files] == 0} { puts "No ZCU104 part0_pins.xml found - search the board_files folder for 'zcu104'"; return }
foreach f $files {
    puts "== $f"
    set fh [open $f]; set txt [read $fh]; close $fh
    foreach line [split $txt "\n"] {
        if {[regexp -nocase {led|dip|sw|clk_125|cpu_reset} $line] && [regexp -nocase {loc=} $line]} { puts [string trim $line] }
    }
}
