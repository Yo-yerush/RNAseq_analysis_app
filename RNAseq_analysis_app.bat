@echo off
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0app_140526"

echo RNAseq dashboard launcher
echo =========================
echo.

rem You can manually set RSCRIPT here if automatic detection fails, for example:
rem set "RSCRIPT=C:\Program Files\R\R-4.5.0\bin\Rscript.exe"
set "RSCRIPT="

for /f "delims=" %%F in ('where Rscript.exe 2^>nul') do (
  if not defined RSCRIPT set "RSCRIPT=%%F"
)

if not defined RSCRIPT (
  if exist "%ProgramFiles%\R" (
    for /f "delims=" %%F in ('where /r "%ProgramFiles%\R" Rscript.exe 2^>nul') do (
      set "RSCRIPT=%%F"
    )
  )
)

if not defined RSCRIPT (
  if exist "%ProgramFiles(x86)%\R" (
    for /f "delims=" %%F in ('where /r "%ProgramFiles(x86)%\R" Rscript.exe 2^>nul') do (
      set "RSCRIPT=%%F"
    )
  )
)

if not defined RSCRIPT (
  if exist "%LOCALAPPDATA%\Programs\R" (
    for /f "delims=" %%F in ('where /r "%LOCALAPPDATA%\Programs\R" Rscript.exe 2^>nul') do (
      set "RSCRIPT=%%F"
    )
  )
)

if not defined RSCRIPT (
  echo ERROR: Rscript.exe was not found.
  echo If R is installed, edit this BAT file and set RSCRIPT manually near the top.
  pause
  exit /b 1
)

echo Found Rscript:
echo   "%RSCRIPT%"
echo.
"%RSCRIPT%" launch_app.R
pause
