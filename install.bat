@echo off
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0app_140526"

echo RNAseq dashboard launcher
echo =========================
echo.

rem You can manually set RSCRIPT here if automatic detection fails, for example:
rem set "RSCRIPT=C:\Program Files\R\R-4.5.0\bin\Rscript.exe"
set "RSCRIPT="

rem 1. Try PATH first
for /f "delims=" %%F in ('where Rscript.exe 2^>nul') do (
  if not defined RSCRIPT set "RSCRIPT=%%F"
)

rem 2. Try common R installation folders under Program Files
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

rem 3. Try user-local installation folder, if present
if not defined RSCRIPT (
  if exist "%LOCALAPPDATA%\Programs\R" (
    for /f "delims=" %%F in ('where /r "%LOCALAPPDATA%\Programs\R" Rscript.exe 2^>nul') do (
      set "RSCRIPT=%%F"
    )
  )
)

if not defined RSCRIPT (
  echo ERROR: Rscript.exe was not found.
  echo.
  echo R may be installed, but this launcher cannot find it.
  echo Please check that this file exists, for example:
  echo   C:\Program Files\R\R-4.x.x\bin\Rscript.exe
  echo.
  echo If it exists, edit this BAT file and set RSCRIPT manually near the top.
  echo Example:
  echo   set "RSCRIPT=C:\Program Files\R\R-4.5.0\bin\Rscript.exe"
  echo.
  pause
  exit /b 1
)

echo Found Rscript:
echo   "%RSCRIPT%"
echo.

echo Checking and installing missing R packages...
"%RSCRIPT%" install_packages.R
if errorlevel 1 (
  echo.
  echo Package installation failed. Check the error above.
  pause
  exit /b 1
)

echo.
echo Launching RNAseq dashboard...
"%RSCRIPT%" launch_app.R
pause
