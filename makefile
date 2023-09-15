# Set compiler, C Language Standard and name of the program
# To use Intel compiler (icx) in windows, insert at terminal before call make:
# call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

CC=clang
CSTD=c17
NAME=CHAOS

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

# Identify flags based on Compiler (icx): 
CFLAGS   :=
ifeq ($(CC),icx)
ifeq ($(OSFLAG), win)
	CFLAGS+=/Qstd:$(CSTD) /Qopenmp /O3 /Qipo /Wall /Werror
else
	CFLAGS+=-std=$(CSTD) -qopenmp -O3 -ipo -Wall -Werror
endif
else ifeq ($(CC),gcc)
	CFLAGS+=-std=$(CSTD) -fopenmp -O3 -Wall -Werror -Wpedantic
else ifeq ($(CC),clang)
	CFLAGS+=-std=$(CSTD) -fopenmp -O3 -Wall -Werror -Wpedantic
else ifeq ($(CC),icc)
	CFLAGS+=-std=$(CSTD) -qopenmp -O3 -Wall -Werror -no-multibyte-chars -diag-disable=10441
endif

# Define Separator, binary extension and execution command based on Operating System
ifeq ($(OSFLAG),win)
	SEP=\\
	EXT=.exe
	EXEC=
else
	SEP=/
	EXT=
	EXEC=./
endif

# Define Important Directories
SRCDIR=src$(SEP)
MODDIR=$(SRCDIR)modules$(SEP)
LIBDIR=$(SRCDIR)libs$(SEP)
BINDIR=bin$(SEP)
OBJDIR=$(BINDIR)obj$(SEP)
OLIBDIR=$(OBJDIR)libs$(SEP)
OMODDIR=$(OBJDIR)modules$(SEP)

# List of .c source files
LIBS=$(wildcard $(LIBDIR)*.c)
MODULES=$(wildcard $(MODDIR)*.c)
SRCS=$(SRCDIR)main.c $(LIBS) $(MODULES) 

# Generate a list of object files from the source files
OBJS=$(patsubst $(SRCDIR)%.c,$(OBJDIR)%.o,$(SRCS))

# The main call of make
all: create_dir $(BINDIR)$(NAME)$(EXT)

# Rule that creates directories based on the Operating System
create_dir:
ifeq ($(OSFLAG),win)
	@if not exist "$(BINDIR)" mkdir "$(BINDIR)"
	@if not exist "$(OBJDIR)" mkdir "$(OBJDIR)"
	@if not exist "$(OLIBDIR)" mkdir "$(OLIBDIR)"
	@if not exist "$(OMODDIR)" mkdir "$(OMODDIR)"
else
	@mkdir -p $(BINDIR)
	@mkdir -p $(OBJDIR)
	@mkdir -p $(OLIBDIR)
	@mkdir -p $(OMODDIR)
endif

# Rule that produce the final binary that depends on the object files
$(BINDIR)$(NAME)$(EXT): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# Rule that produce the object files. Each object file depends on its corresponding source file and the header file
$(OBJDIR)%.o: $(SRCDIR)%.c $(wildcard $(SRCDIR)*.h)
	$(CC) $(CFLAGS) -c $< -o $@

run: $(BINDIR)$(NAME)$(EXT)
	@$(EXEC)$^

plot:
ifeq ($(OSFLAG), win)
	@scripts\plot.bat
else
	@bash scripts/plot.sh
endif

test:
ifeq ($(OSFLAG), win)
	@$(CC) -o tests\\test tests\\test.c
	@tests\\test.exe
else
	@$(CC) -o tests/test tests/test.c
	@./tests/test
endif

clean: 
ifeq ($(OSFLAG), win)
	@del $(OBJDIR)*.o
	@del $(OLIBDIR)*.o
	@del $(OMODDIR)*.o
	@del $(BINDIR)$(NAME)$(EXT)
else
	@rm -rf $(OBJDIR)*.o
	@rm -rf $(OLIBDIR)*.o
	@rm -rf $(OMODDIR)*.o
	@rm -f $(BINDIR)$(NAME)$(EXT)
endif
