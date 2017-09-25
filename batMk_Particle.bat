@echo off
setlocal enabledelayedexpansion

REM  put on .exe files. 
REM  "../data/params.txt" must exist upon .exe files' directory 

for /f "skip=1" %%a in (..\data\params.txt) do (
  call %1 %%a
)
