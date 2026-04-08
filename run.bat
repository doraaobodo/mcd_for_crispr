
@echo off
setlocal enabledelayedexpansion

:: ------------------------------------------
:: Configuration
:: ------------------------------------------
set "SCRIPT_DIR=%~dp0"
set "R_SCRIPT_NAME=main.R"
set "MAIN_R=%SCRIPT_DIR%%R_SCRIPT_NAME%"
set "RSCRIPT_EXE="

:: ------------------------------------------
:: Check that the R script exists
:: ------------------------------------------
if not exist "%MAIN_R%" (
    echo Error: Could not find %R_SCRIPT_NAME% in:
    echo %SCRIPT_DIR%
    echo.
    pause
    exit /b 1
)

:: ------------------------------------------
:: First try to find Rscript on PATH
:: ------------------------------------------
where Rscript >nul 2>nul
if not errorlevel 1 (
    for /f "delims=" %%I in ('where Rscript') do (
        set "RSCRIPT_EXE=%%I"
        goto :run
    )
)

:: ------------------------------------------
:: Otherwise search common install locations
:: Newest version first
:: ------------------------------------------

call :find_latest_rscript "C:\Program Files\R"
if defined RSCRIPT_EXE goto :run

call :find_latest_rscript "C:\Program Files (x86)\R"
if defined RSCRIPT_EXE goto :run

:: ------------------------------------------
:: If still not found, fail clearly
:: ------------------------------------------
echo Error: Rscript was not found on this system.
echo Please install R and make sure Rscript is available.
echo See the README.md for setup instructions.
echo.
pause
exit /b 1

:: ------------------------------------------
:: Run the R script
:: ------------------------------------------
:run
echo ==========================================
echo MCD UI Launcher
echo ==========================================
echo Using R at:
echo %RSCRIPT_EXE%
echo.
echo Running script:
echo %MAIN_R%
echo ==========================================
echo.

"%RSCRIPT_EXE%" --vanilla "%MAIN_R%" %*
set "EXITCODE=%ERRORLEVEL%"

echo.
if errorlevel 1 (
    echo The script ended with an error.
) else (
    echo The script finished successfully.
)
pause
exit /b %EXITCODE%

:: ------------------------------------------
:: Function: find latest Rscript in a root dir
:: Looks for newest R-* folder first
:: ------------------------------------------
:find_latest_rscript
set "R_ROOT=%~1"

if not exist "%R_ROOT%" goto :eof

for /f "delims=" %%I in ('dir /b /ad /o-n "%R_ROOT%\R-*" 2^>nul') do (
    if exist "%R_ROOT%\%%I\bin\Rscript.exe" (
        set "RSCRIPT_EXE=%R_ROOT%\%%I\bin\Rscript.exe"
        goto :eof
    )
    if exist "%R_ROOT%\%%I\bin\x64\Rscript.exe" (
        set "RSCRIPT_EXE=%R_ROOT%\%%I\bin\x64\Rscript.exe"
        goto :eof
    )
)

goto :eof