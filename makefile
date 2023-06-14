# Set compiler, C Language Standard and name of the program
# To use Intel compiler (icx) in windows, insert at terminal before call make:
# call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

CC=gcc
CSTD=c17
NAME=CHAOS
CFORGENAME = CHAOSForge
# Indentify Operating System
OSFLAG 				:=
ifeq ($(OS),Windows_NT)
	OSFLAG+=win
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OSFLAG+=linux
	endif
	ifeq ($(UNAME_S),Darwin)
		OSFLAG+=macos
	endif
endif
# Identify flags based on Compiler (icx): /Wall /Werror     (gcc): -Wall -Werror -Wpedantic
CFLAGS   :=
ifeq ($(CC),icx) 
ifeq ($(OSFLAG), win)
	CFLAGS+=/Qstd:$(CSTD) /Qopenmp /O3 /Qipo 
else
	CFLAGS+=-std=$(CSTD) -qopenmp -O3 -ipo 
endif
else ifeq ($(CC),gcc)
	CFLAGS+=-std=$(CSTD) -fopenmp -O3 
else ifeq ($(CC),clang)
	CFLAGS+=-std=$(CSTD) -fopenmp -O3
else ifeq ($(CC),icc)
	CFLAGS+=-std=$(CSTD) -qopenmp -no-multibyte-chars -O3 -ipo -diag-disable=10441
endif

# -Wall -Werror -Wpedantic

#define .c files to be compiled
LIBS=src/libs/odesystems.c src/libs/interface.c src/libs/iofiles.c src/libs/nldyn.c src/libs/nlosc.c src/libs/customcalc.c src/libs/basic.c src/libs/odesolvers.c src/libs/msg.c
MODULES=src/modules/convergence_test.c src/modules/time_series.c src/modules/poinc_map.c src/modules/lyap_exp_wolf.c src/modules/ftime_series.c src/modules/bifurcation.c src/modules/fbifurcation.c src/modules/dyndiag.c src/modules/fdyndiag.c src/modules/epbasin.c src/modules/forcedbasin.c 
HOSMODULES=src/modules/HOS_time_series.c src/modules/HOS_ftime_series.c src/modules/HOS_bifurcation.c src/modules/HOS_fbifurcation.c src/modules/HOS_dyndiag.c src/modules/HOS_fdyndiag.c src/modules/HOS_fforcedbasin.c
FILES=src/main.c $(LIBS) $(MODULES) $(HOSMODULES) 

CFORGEFILES = src/libs/msg.c src/libs/iofiles.c src/libs/basic.c src/cforge.c
# Identify CFORGE flags based on Compiler /Wall /Werror
CFORGEFLAGS   :=
ifeq ($(CC),icx) 
ifeq ($(OSFLAG), win)
	CFORGEFLAGS+=/Qstd:$(CSTD) 
else
	CFORGEFLAGS+=-std=$(CSTD) 
endif
else ifeq ($(CC),gcc)
	CFORGEFLAGS+=-std=$(CSTD)
else ifeq ($(CC),clang)
	CFORGEFLAGS+=-std=$(CSTD)
else ifeq ($(CC),icc)
	CFORGEFLAGS+=-std=$(CSTD) -no-multibyte-chars -diag-disable=10441
endif

all:
	@echo $(OS)
ifeq ($(OSFLAG), win)
	@if not exist "bin\" mkdir "bin\"
	$(CC) $(CFLAGS) -o bin\$(NAME) $(FILES)
#	@scripts\assign_icon_win.bat $(NAME)
else ifeq ($(OSFLAG), macos)
	@mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$(NAME) $(FILES)
#	@bash scripts/assign_icon_macos.sh
else
	@mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$(NAME) $(FILES)
endif

run: 
ifeq ($(OSFLAG), win)
	@bin\$(NAME).exe
else
	@./bin/$(NAME)
endif

plot:
ifeq ($(OSFLAG), win)
	@scripts\plot.bat
else
	@bash scripts/plot.sh
endif

cforge:
	@echo $(OS)
ifeq ($(OSFLAG), win)
	@if not exist "bin\" mkdir "bin\"
	$(CC) $(CFORGEFLAGS) -o bin\$(CFORGENAME) $(CFORGEFILES)
#	@scripts\assign_icon_win.bat $(NAME)
else ifeq ($(OSFLAG), macos)
	@mkdir -p bin
	$(CC) $(CFORGEFLAGS) -o bin/$(CFORGENAME) $(CFORGEFILES)
#	@bash scripts/assign_icon_macos.sh
else
	@mkdir -p bin
	$(CC) $(CFORGEFLAGS) -o bin/$(CFORGENAME) $(CFORGEFILES)
endif

run_cforge:
ifeq ($(OSFLAG), win)
	@bin\$(CFORGENAME).exe
else
	@./bin/$(CFORGENAME)
endif

test_cforge: cforge run_cforge

RIFLAGS   :=
ifeq ($(CC),icx) 
ifeq ($(OSFLAG), win)
	RIFLAGS+=/Qstd:$(CSTD) /Wall /Werror
else
	RIFLAGS+=-std=$(CSTD) 
endif
else ifeq ($(CC),gcc)
	RIFLAGS+=-std=$(CSTD)
else ifeq ($(CC),clang)
	RIFLAGS+=-std=$(CSTD)
else ifeq ($(CC),icc)
	RIFLAGS+=-std=$(CSTD) -no-multibyte-chars -diag-disable=10441
endif

RINAME=read_in
RIFILES=src/libs/msg.c src/libs/iofiles.c src/libs/basic.c src/read_input.c
ri:
ifeq ($(OSFLAG), win)
	@if not exist "bin\" mkdir "bin\"
	$(CC) $(RIFLAGS) -o bin\$(RINAME) $(RIFILES)
#	@scripts\assign_icon_win.bat $(NAME)
else ifeq ($(OSFLAG), macos)
	@mkdir -p bin
	$(CC) $(RIFLAGS) -o bin/$(RINAME) $(RIFILES)
#	@bash scripts/assign_icon_macos.sh
else
	@mkdir -p bin
	$(CC) $(RIFLAGS) -o bin/$(RINAME) $(RIFILES)
endif

run_ri:
ifeq ($(OSFLAG), win)
	@bin\$(RINAME).exe
else
	@./bin/$(RINAME)
endif

test_ri: ri run_ri

clean_bin:
ifeq ($(OSFLAG), win)
	@del /Q bin
else
	@rm -r bin
endif

clean: 
ifeq ($(OSFLAG), win)
	@del /Q bin
else
	@rm -r bin
endif