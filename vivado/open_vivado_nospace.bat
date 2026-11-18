@echo off
REM Test launcher: starts Vivado with its settings folder on a path WITHOUT spaces (the Windows user name has spaces).
REM Close Vivado completely first (and end any hw_server / xsdb / cs_server processes in Task Manager).
if not exist C:\vivado_cfg\AppData mkdir C:\vivado_cfg\AppData
set APPDATA=C:\vivado_cfg\AppData
call C:\AMDDesignTools\2025.2\Vivado\settings64.bat
start "" vivado C:\SpaceGuard\SpaceGuard_FPGA-Based-Space-Radiation-Monitoring-and-Anomaly-Detection\rtl\rtl\rtl.xpr
