@echo off
setlocal

REM Save original directory
set "ORIG_DIR=%cd%"

REM Change to the license server tools directory
cd "C:\Program Files (x86)\Accelrys\LicensePack\win32\bin" 

REM Define temp file
set "OUTFILE=%TEMP%\ms_visualizer_usage.txt"

REM Clear previous output
> "%OUTFILE%" echo Materials Studio Visualizer License Usage
>>"%OUTFILE%" echo Timestamp: %date% %time%
>>"%OUTFILE%" echo.

REM Run lmutil and capture MS_visualizer section
for /f "tokens=*" %%A in ('lmutil lmstat -a') do (
    echo %%A | findstr /C:"Users of MS_visualizer" >nul
    if not errorlevel 1 (
        >>"%OUTFILE%" echo %%A
        set "capture=1"
    ) else (
        echo %%A | findstr /C:"Users of " >nul
        if not errorlevel 1 (
            set "capture="
        )
    )

    if defined capture (
        >>"%OUTFILE%" echo %%A
    )
)

REM Return to original directory
cd "%ORIG_DIR%"

REM Open the temp file in Notepad
notepad "%OUTFILE%"

endlocal
