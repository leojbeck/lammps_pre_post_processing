@echo off
REM Batch script to display current users of MS_visualizer from FlexLM

REM Run lmutil and capture output
cd "C:\Program Files (x86)\Accelrys\LicensePack\win32\bin" 
for /f "tokens=*" %%A in ('lmutil lmstat -a') do (
    echo %%A | findstr /C:"Users of MS_visualizer" >nul
    if not errorlevel 1 (
        echo %%A
        set "capture=1"
    ) else (
        echo %%A | findstr /C:"Users of " >nul
        if not errorlevel 1 (
            set "capture="
        )
    )

    if defined capture (
        echo %%A
    )
)

echo.
echo Done.
pause
