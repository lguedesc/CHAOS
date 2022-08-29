call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

::icx /Qstd:c17 /Zi /DEBUG /Qopenmp /O3 /Ot /Qipo src/main.c src/libs/edosystems.c src/libs/nldyn.c src/libs/iofiles.c src/libs/interface.c src/modules/time_series.c src/modules/poinc_map.c src/modules/lyap_exp_wolf.c src/modules/ftime_series.c src/modules/bifurcation.c src/modules/fbifurcation.c src/modules/dyndiag.c src/modules/epbasin.c src/modules/forcedbasin.c
icx /Qstd:c17 /MD /Qopenmp /O3 /Ot /Qipo src/main.c src/libs/edosystems.c src/libs/nldyn.c src/libs/iofiles.c src/libs/interface.c src/modules/time_series.c src/modules/poinc_map.c src/modules/lyap_exp_wolf.c src/modules/ftime_series.c src/modules/bifurcation.c src/modules/fbifurcation.c src/modules/dyndiag.c src/modules/epbasin.c src/modules/forcedbasin.c

@ main.exe
::@ %USERPROFILE%/anaconda3/python.exe -B -m TimeSeries.out.plot_rk4
