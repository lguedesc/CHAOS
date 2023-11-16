#ifndef MSG_H
#define MSG_H

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

#define RESET_STYLE        "\x1b[0m"

#define COLOR_BLACK        "\x1b[30m"
#define COLOR_RED          "\x1b[31m"
#define COLOR_GREEN        "\x1b[32m"
#define COLOR_YELLOW       "\x1b[33m"
#define COLOR_BLUE         "\x1b[34m"
#define COLOR_MAGENTA      "\x1b[35m"
#define COLOR_CYAN         "\x1b[36m"
#define COLOR_WHITE        "\x1b[37m"

#define BACKGROUND_BLACK   "\x1b[40m"
#define BACKGROUND_RED     "\x1b[41m"
#define BACKGROUND_GREEN   "\x1b[42m"
#define BACKGROUND_YELLOW  "\x1b[43m"
#define BACKGROUND_BLUE    "\x1b[44m"
#define BACKGROUND_MAGENTA "\x1b[45m"
#define BACKGROUND_CYAN    "\x1b[46m"
#define BACKGROUND_WHITE   "\x1b[47m"

#define STYLE_BOLD         "\x1b[1m"
#define STYLE_ITALIC       "\x1b[3m"
#define STYLE_UNDERLINE    "\x1b[4m"

/* Colors */
void yellow(void);
void green(void);
void red(void);
void blue(void);
void magenta(void);
void cyan(void);
void reset_color(void);

/* Generic Messagens */
void print_warning(char *msg, ...);
void print_success(char *msg, ...);
void print_error(char *msg, ...);
void print_blue(char *msg, ...);
void print_magenta(char *msg, ...);
void print_cyan(char *msg, ...);
void print_debug(char *msg, ...);
/* Specific Messages */
void print_exit_prog(void);
void reset_program(char *issue);

#endif
