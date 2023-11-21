echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v232\fluent/ntbin/win64/winkill.exe"

start "tell.exe" /B "C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v232\fluent\ntbin\win64\tell.exe" LAPTOP-57O9LBH0 19281 CLEANUP_EXITING
timeout /t 1
"C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v232\fluent\ntbin\win64\kill.exe" tell.exe
if /i "%LOCALHOST%"=="LAPTOP-57O9LBH0" (%KILL_CMD% 4504) 
if /i "%LOCALHOST%"=="LAPTOP-57O9LBH0" (%KILL_CMD% 8468) 
if /i "%LOCALHOST%"=="LAPTOP-57O9LBH0" (%KILL_CMD% 12816) 
if /i "%LOCALHOST%"=="LAPTOP-57O9LBH0" (%KILL_CMD% 12148) 
if /i "%LOCALHOST%"=="LAPTOP-57O9LBH0" (%KILL_CMD% 17748) 
if /i "%LOCALHOST%"=="LAPTOP-57O9LBH0" (%KILL_CMD% 17160)
del "D:\NITCMECHSEM7\Computational Fluid Dynamics\CFD assignment\cleanup-fluent-LAPTOP-57O9LBH0-17748.bat"
