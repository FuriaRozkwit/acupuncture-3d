@echo off
chcp 65001 >nul
title Atlas Akupunktury 3D - Konsola Kliniczna
cd /d "%~dp0"
echo ========================================================
echo   ATLAS AKUPUNKTURY 3D - POSAZEK Z GLINY I SIATKA CUNOWA
echo ========================================================
echo   [1/2] Inicjalizacja silnika 3D i siatki anatomicznej...
echo   [2/2] Okno graficzne pojawi sie za chwile.
echo.
echo   UWAGA: Nie zamykaj tego okna terminala!
echo   Po kliknieciu punktu w oknie 3D tutaj pojawia sie pelny
echo   opis kliniczny, 8 Naczyn Nadzwyczajnych i pary Yuan-Luo.
echo ========================================================
echo.
set DISTRO=Ubuntu-24.04
wsl.exe -d %DISTRO% -u root bash tools/run_wsl.sh --surface %*
if errorlevel 1 (
    echo.
    echo [BLAD] Wystapil problem z uruchomieniem okna graficznego.
    pause
)
