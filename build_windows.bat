@echo off
setlocal

REM One-click Windows build script for RNA Forge 5' Capping Efficiency Tool.
REM Run this from the project folder after installing requirements:
REM     python -m pip install -r requirements.txt

python -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onefile ^
  --windowed ^
  --name "RNA_Forge_Capping_Efficiency_Tool" ^
  main.py

if errorlevel 1 (
  echo Build failed.
  exit /b 1
)

echo Build complete. EXE should be in the dist folder.
endlocal
