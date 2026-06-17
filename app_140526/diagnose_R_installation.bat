@echo off
setlocal EnableExtensions EnableDelayedExpansion

echo Checking R installation
echo =======================
echo.

echo 1. Searching PATH:
where R.exe
where Rscript.exe

echo.
echo 2. Searching Program Files R folders:
if exist "%ProgramFiles%\R" (
  where /r "%ProgramFiles%\R" Rscript.exe
) else (
  echo No folder: "%ProgramFiles%\R"
)

echo.
if exist "%ProgramFiles(x86)%\R" (
  where /r "%ProgramFiles(x86)%\R" Rscript.exe
) else (
  echo No folder: "%ProgramFiles(x86)%\R"
)

echo.
echo 3. Searching user-local R folder:
if exist "%LOCALAPPDATA%\Programs\R" (
  where /r "%LOCALAPPDATA%\Programs\R" Rscript.exe
) else (
  echo No folder: "%LOCALAPPDATA%\Programs\R"
)

echo.
pause
