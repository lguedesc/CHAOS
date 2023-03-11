CCwin = icx
CCunix = icc
winCFLAGS = /Qstd:c17 /MD /Qopenmp /O3 /Ot /Qipo

unixCFLAGS = -std=c17 -no-multibyte-chars -qopenmp -O3

FILES = src/main.c src/libs/odesystems.c src/libs/interface.c src/libs/iofiles.c src/libs/nldyn.c src/libs/nlosc.c src/libs/customcalc.c src/modules/convergence_test.c src/modules/time_series.c src/modules/poinc_map.c src/modules/lyap_exp_wolf.c src/modules/ftime_series.c src/modules/bifurcation.c src/modules/fbifurcation.c src/modules/dyndiag.c src/modules/fdyndiag.c src/modules/epbasin.c src/modules/forcedbasin.c src/modules/OS_time_series.c src/modules/OS_ftime_series.c src/modules/OS_bifurcation.c src/modules/OS_fbifurcation.c src/modules/OS_dyndiag.c src/modules/OS_fdyndiag.c src/modules/OS_fforcedbasin.c
NAME = CHAOS

initwin:
	@call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

windows:
	@$(CCwin) $(winCFLAGS) $(FILES) -o bin/$(NAME)
	@bin\CHAOS.exe

unix:
	@$(CCunix) $(unixCFLAGS) $(FILES) -o bin/$(NAME)
	@./bin/CHAOS

test:
	@$(CC) $(CFLAGS) test/test.c -o tests
	@./tests

clean_bin:
	@rm -f bin/*

clean:
	@rm -f bin/* out/*