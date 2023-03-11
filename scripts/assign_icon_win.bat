@echo off
:: Create program.rc file
echo IDI_ICON_1 ICON "assets\icons\chaos.ico" > program.rc
:: Generate the resource file
rc program.rc
:: Create bin folder if it doesnt exist
if not exist "bin\" mkdir bin\
:: Link the created object with the created resource file
link bin\%1.obj program.res -out:bin/%1.exe
:: Remove non-wanted files
del bin\%1.obj
del program.res
del program.rc