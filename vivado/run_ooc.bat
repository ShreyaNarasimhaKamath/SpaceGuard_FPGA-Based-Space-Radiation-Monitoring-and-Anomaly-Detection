@echo off
rem Run from anywhere: runs the SpaceGuard out-of-context Vivado flow. Edit PART for another ZCU board.
set PART=%1
if "%PART%"=="" set PART=xczu7ev-ffvc1156-2-e
call C:\AMDDesignTools\2025.2\Vivado\settings64.bat
cd /d "%~dp0.."
vivado -mode batch -source vivado\run_ooc.tcl -tclargs %PART% -log results\vivado\vivado_%PART%.log -journal results\vivado\vivado_%PART%.jou
pause
