@echo off

call "E:\Visual Studio 2017\Common7\Tools\..\..\VC\Auxiliary\Build\VCVARSALL.BAT" AMD64


call  "\\Yudi\E$\MatLab R2017b\bin\win64\checkMATLABRootForDriveMap.exe" "\\Yudi\E$\MatLab R2017b"  > mlEnv.txt
for /f %%a in (mlEnv.txt) do set "%%a"\n
cd .

if "%1"=="" (nmake MATLAB_ROOT=%MATLAB_ROOT% ALT_MATLAB_ROOT=%ALT_MATLAB_ROOT% MATLAB_BIN=%MATLAB_BIN% ALT_MATLAB_BIN=%ALT_MATLAB_BIN%  -f buttordBandpassOnly_rtw.mk all) else (nmake MATLAB_ROOT=%MATLAB_ROOT% ALT_MATLAB_ROOT=%ALT_MATLAB_ROOT% MATLAB_BIN=%MATLAB_BIN% ALT_MATLAB_BIN=%ALT_MATLAB_BIN%  -f buttordBandpassOnly_rtw.mk %1)
@if errorlevel 1 goto error_exit

exit /B 0

:error_exit
echo The make command returned an error of %errorlevel%
An_error_occurred_during_the_call_to_make
