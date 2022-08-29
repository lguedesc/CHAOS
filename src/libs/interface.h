#include <stdio.h>
#include <stdlib.h>

void clear_screen();
void partition(int mode, size_t maxlength);
void welcome_header(size_t maxlength);
void option_header(char *tbnames, size_t maxlength);
void mixed_header(char *names1, char *names2, size_t maxlength);
unsigned int choose_option(char **options, size_t n, size_t maxlength, char *category);
void invalid_option(unsigned int option, char* category, size_t maxlength);
void end_of_execution(size_t maxlength);
void identify_simulation(unsigned int toolbox, unsigned int *system, unsigned int *module, char **toolboxesNames, char **systemNames, char** moduleNames, size_t numofsystems, size_t maxlength, size_t numofmodules);