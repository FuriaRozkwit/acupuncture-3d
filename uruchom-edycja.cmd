@echo off
chcp 65001 >nul
title Atlas Akupunktury 3D - Tryb Edycji i Profile
cd /d "%~dp0"
echo ========================================================
echo   ATLAS AKUPUNKTURY 3D - TRYB EDYCJI I PROFILE
echo ========================================================
echo   [1/2] Inicjalizacja silnika 3D...
echo   [2/2] Okno pojawi sie za chwile.
echo.
echo   Skroty w oknie:
echo     [t] - Wlacz / wylacz tryb edycji (kliknij punkt by przesuwac)
echo     [v] - Menu wyboru profili i konfiguracji
echo ========================================================
echo.
set DISTRO=Ubuntu-24.04
wsl.exe -d %DISTRO% -u root bash tools/run_edit_wsl.sh %*
if errorlevel 1 (
    echo.
    echo [BLAD] Wystapil problem z uruchomieniem okna edytora.
    pause
)
