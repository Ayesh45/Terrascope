@echo off
setlocal
if "%EE_PROJECT%"=="" (
  echo Please set EE_PROJECT to your Google Earth Engine-enabled Cloud project ID first.
  echo Example in this terminal: set EE_PROJECT=your-project-id
  pause
  exit /b 1
)
if not exist .venv\Scripts\python.exe (
  echo Creating Python virtual environment...
  py -3.11 -m venv .venv
  if errorlevel 1 py -m venv .venv
  .venv\Scripts\python.exe -m pip install --upgrade pip
  .venv\Scripts\python.exe -m pip install -r requirements.txt
)
.venv\Scripts\python.exe -c "import ee; ee.Authenticate()"
.venv\Scripts\python.exe app.py
