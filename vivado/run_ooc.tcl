# SpaceGuard: out-of-context synthesis + implementation of sg_top variants on a Zynq UltraScale+ part.
# Usage A (Vivado GUI Tcl console):  source C:/SpaceGuard/SpaceGuard_FPGA-Based-Space-Radiation-Monitoring-and-Anomaly-Detection/vivado/run_ooc.tcl
# Usage B (batch):  vivado -mode batch -source vivado/run_ooc.tcl -tclargs <part>
# Default part: xczu7ev-ffvc1156-2-e (ZCU104/106). Others (need a license for the big parts): xczu9eg-ffvb1156-2-e (ZCU102), xczu28dr-ffvg1517-2-e (ZCU111).
set part "xczu7ev-ffvc1156-2-e"
if {[info exists ::argv] && [llength $::argv] > 0} { set part [lindex $::argv 0] }
# GUI: an open project (e.g. rtl.xpr) blocks create_project -in_memory, so remember it, close it, and reopen it at the end.
set prev_xpr ""
if {![catch {current_project} cp] && $cp ne ""} {
  set prev_xpr [file join [get_property DIRECTORY [current_project]] ${cp}.xpr]
  close_project
}
set root [file normalize [file join [file dirname [info script]] ..]]
set outd [file join $root results vivado]; file mkdir $outd
set srcs [list rtl/anomaly_detector/sg_next.v rtl/anomaly_detector/sg_detector.v rtl/anomaly_detector/sg_cdet.v rtl/top/sg_ctrl.v rtl/top/sg_top.v]
# name  KIND DET_TMR CTRL_TMR
set variants {
  {bare_ema   0 0 0}
  {tmr_ema    0 1 1}
  {bare_count 1 0 0}
  {tmr_count  1 1 1}
}
set csv [open [file join $outd "summary_$part.csv"] w]
puts $csv "part,variant,LUT,FF,CARRY8,BRAM,DSP,WNS_ns,period_ns,fmax_MHz_est,total_power_W,dynamic_W,static_W"
foreach v $variants {
  lassign $v name kind dtmr ctmr
  puts "=== $name ($part) ==="
  create_project -in_memory -part $part
  foreach f $srcs { read_verilog [file join $root $f] }
  read_xdc -mode out_of_context [file join $root vivado ooc.xdc]
  synth_design -top sg_top -part $part -mode out_of_context -flatten_hierarchy none \
      -generic N=2 -generic K=2 -generic KIND=$kind -generic DET_TMR=$dtmr -generic CTRL_TMR=$ctmr
  opt_design; place_design; route_design
  set pre [file join $outd "${name}_$part"]
  report_utilization -file ${pre}_util.rpt
  report_timing_summary -file ${pre}_timing.rpt
  report_power -file ${pre}_power.rpt
  set util [report_utilization -return_string]
  proc grab {txt pat} { if {[regexp $pat $txt -> val]} {return $val}; return NA }
  set lut  [grab $util {\| CLB LUTs\s*\|\s*(\d+)}]
  set ff   [grab $util {\| CLB Registers\s*\|\s*(\d+)}]
  set car  [grab $util {\| CARRY8\s*\|\s*(\d+)}]
  set bram [grab $util {\| Block RAM Tile\s*\|\s*(\d+)}]
  set dsp  [grab $util {\| DSPs\s*\|\s*(\d+)}]
  set wns  [get_property SLACK [get_timing_paths -max_paths 1 -setup]]
  set per 3.333
  set fmax [expr {$wns eq "" ? "NA" : 1000.0/($per - $wns)}]
  set pw [report_power -return_string]
  set tp [grab $pw {Total On-Chip Power \(W\)\s*\|\s*([0-9.]+)}]
  set dp [grab $pw {Dynamic \(W\)\s*\|\s*([0-9.]+)}]
  set sp [grab $pw {Device Static \(W\)\s*\|\s*([0-9.]+)}]
  puts $csv "$part,$name,$lut,$ff,$car,$bram,$dsp,$wns,$per,$fmax,$tp,$dp,$sp"; flush $csv
  # how many state flops survived (TMR copies must not be merged)
  set n [llength [get_cells -hier -filter {IS_SEQUENTIAL && REF_NAME =~ FD*}]]
  puts "flops after impl: $n"
  close_project
}
close $csv
puts "DONE -> $outd"
if {$prev_xpr ne "" && [file exists $prev_xpr]} { open_project $prev_xpr }
