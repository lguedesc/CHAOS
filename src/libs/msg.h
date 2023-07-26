#ifndef MSG_H
#define MSG_H

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

/* Generic Messagens */
void print_warning(char *msg, ...);
void print_success(char *msg, ...);
void print_error(char *msg, ...);
void print_blue(char *msg, ...);
void print_purple(char *msg, ...);
void print_cyan(char *msg, ...);
void print_debug(char *msg, ...);
/* Specific Messages */
void print_exit_prog(void);
void reset_program(char *issue);

#endif
