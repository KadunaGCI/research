@echo off
setlocal enabledelayedexpansion

REM  put on .exe files. 
REM  "../data" and "../result" must exist upon .exe files' directory 

for /f "skip=1" %%a in (..\data\params.txt) do (

rem get time
set tmp_date=!date!
set tmp_time=!time!
echo !tmp_time!
set year=!tmp_date:~0,4!
set month=!tmp_date:~5,2!
set day=!tmp_date:~8,2!
REM kunikunosaku
set hour=!tmp_time:~1,1!
set minute=!tmp_time:~3,2!
set umin=!tmp_time:~6,2!
set now=!month!!day!!year!_!hour!!minute!!umin!
set path=..\result\!now!
echo !path!
  
md !path!
call %1 %%a
  
for /f "delims=," %%a in ('dir ..\data\output* /b /s') do (
REM  echo %%a
  copy "%%a" !path!
)
for /f "delims=," %%a in ('dir ..\data\particle* /b /s') do (
  copy "%%a" !path!
)
echo --------------------------------------------------------------------
)

pause
