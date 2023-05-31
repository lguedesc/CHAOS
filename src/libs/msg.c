#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

/* Colors */
void yellow() {
    printf("\033[0;33m");
}
void green() {
    printf("\033[0;32m");
}
void red() {
    printf("\033[0;31m");
}

void blue() {
    printf("\033[0;34m");
}

void purple() {
    printf("\033[0;35m");
}

void cyan() {
    printf("\033[0;36m");
}
void reset_color() {
    printf("\033[0m");
}

/* Generic Messagens */

void print_warning(char *msg, ...) {
    yellow();
    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
    reset_color();
}
void print_success(char *msg, ...) {
    green();
    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
    reset_color();
}
void print_error(char *msg, ...) {
    red();
    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
    reset_color();
}
void print_blue(char *msg, ...) {
    blue();
    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
    reset_color();
}
void print_purple(char *msg, ...) {
    purple();
    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
    reset_color();
}
void print_cyan(char *msg, ...) {
    cyan();
    va_list args;
    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);
    reset_color();
}

/* Specific Messages */

void print_exit_prog() {
    printf("Exiting program...\n");
}
void reset_program(char *issue) {
    print_error("%s", issue);
    print_error("To solve the issue, try to do the following procedure:\n");
    print_error("   1. Download the CHAOS package again at 'https://github.com/lguedesc/CHAOS.git'\n");
    print_error("   2. Compile and run the program as indicated in the documentation.\n");
    print_error("   3. Compile and run AddSys program as indicated in the documentation.\n");
    print_error("   4. If the problem persists, please open an issue at 'https://github.com/lguedesc/CHAOS/issues' describing the problem.\n");
    print_exit_prog();
    exit(EXIT_FAILURE);
}