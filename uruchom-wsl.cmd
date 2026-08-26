@echo off
rem Uruchamia atlas w WSL2 (okno przez WSLg).
rem Potrzebne, bo Smart App Control na Windows blokuje niepodpisane biblioteki VTK.
rem Pierwsze uruchomienie wymaga: wsl -d Ubuntu-24.04 -u root bash tools/setup_wsl.sh
setlocal
set DISTRO=Ubuntu-24.04
for /f "usebackq delims=" %%p in (`wsl -d %DISTRO% wslpath -a "%~dp0tools/run_wsl.sh"`) do set SCRIPT=%%p
if "%SCRIPT%"=="" (
    echo Nie udalo sie ustalic sciezki do tools/run_wsl.sh w WSL.
    exit /b 2
)
wsl -d %DISTRO% -u root bash "%SCRIPT%" %*
