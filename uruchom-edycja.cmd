@echo off
rem Atlas z edycja punktow: [t] tryb edycji, [v] konfiguracje (fabryczna / profile).
rem Domyslnie laduje sie profil wybrany jako domyslny; fabryczna: uruchom-edycja.cmd --fabryczna
cd /d "%~dp0"
set DISTRO=Ubuntu-24.04
wsl.exe -d %DISTRO% -u root bash tools/run_edit_wsl.sh %*
