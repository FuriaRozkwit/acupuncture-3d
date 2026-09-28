@echo off
rem Uruchamia atlas w WSL2 (okno przez WSLg).
rem Potrzebne, bo Smart App Control na Windows blokuje niepodpisane biblioteki VTK.
cd /d "%~dp0"
set DISTRO=Ubuntu-24.04
wsl.exe -d %DISTRO% -u root bash tools/run_wsl.sh --surface %*

