#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void clear_screen() {
    #ifdef _WIN32
        system("cls");
    #else
        system("clear");
    #endif
}

void partition(int mode, size_t maxlength) {
    if (mode == 1) {
        printf("%c", 218);
        for (size_t i = 0; i < maxlength-2; i++) {
            printf("%c", 196);
        }
        printf("%c", 191);
    }
    else if (mode == 2) {
        printf("%c", 192);
        for (size_t i = 0; i < maxlength-2; i++) {
            printf("%c", 196);
        }
        printf("%c", 217);
    }
    else if (mode == 3) {
        printf("%c",195);
        for (size_t i = 0; i < maxlength-2; i++) {
            printf("%c", 196);
        }
        printf("%c",180);
    }
    else if (mode == 4) {
        printf(" ");
        for (size_t i = 0; i < maxlength-2; i++) {
            printf("%c", 196);
        }
        printf(" ");
    }
    else {
        printf("DEBUG WARNING: Please select 1 for main partition or 2 for secondary partition for the interface programming");
        exit(0);
    }
    printf("\n");
}

void welcome_header(size_t maxlength) {
    //clear_screen();
    char *message = "Welcome to CHAOS, a Nonlinear Dynamics Package for Harmonically Forced Systems";
    partition(1, maxlength);
    int padlen = (maxlength - strlen(message)) / 2;
    printf("%c%*s%s%*c\n", 179, padlen-1, "", message, padlen, 179);
    partition(2, maxlength);
}

void option_header(char *names, size_t maxlength) {
    int padlen = (maxlength - strlen(names)) / 2;
    partition(1, maxlength);
    printf("%*s%s%*s\n", padlen, "", names, padlen, "");
    partition(3, maxlength);
}

unsigned int choose_option(char **options, size_t n, size_t maxlength, char *category) {
    unsigned int choice;
    for (int i = 0; i < n; i++) {
        printf("%s%-3d%s%s\n", "  ", i+1, "- ", options[i]);
    }
    printf("%s%-3d%s\n", "  ", 0, "- EXIT");
    partition(2, maxlength);
    printf("%s%s%s","  Choose a ", category, " for the Analysis: ");
    //scanf("%u", &choice);
    if (scanf(" %u",&choice) == 1){
        return choice;
    }
    else {
        return choice;
    }
}

void mixed_header(char *names1, char *names2, size_t maxlength) {
    int padlen = (maxlength - strlen(names1) - strlen(names2)) / 2;
    printf("%*s%s%s%s%*s\n", padlen, "", names1, " - ", names2 , padlen, "");
    partition(2, maxlength);
}

void invalid_option(unsigned int option, char* category, size_t maxlength) {
    printf("%s%s%s%u%s", "  ", category, " ", option, " : Invalid option. Press any key to exit program...\n");
    partition(2, maxlength);
    while(getchar()!='\n'){}
    getchar(); // wait for ENTER
    exit(0);
}

void end_of_execution(size_t maxlength) {
    printf("%s", "\n  Execution ended successfully! Press any key to exit program...\n");
    partition(4, maxlength);
    while(getchar()!='\n'){}
    getchar(); // wait for ENTER
}

void identify_simulation(unsigned int toolbox, unsigned int *system, unsigned int *module, char **toolboxesNames, char **systemNames, char** moduleNames, size_t numofsystems, size_t maxlength, size_t numofmodules) {
    clear_screen();
    welcome_header(maxlength);
    option_header(toolboxesNames[toolbox-1], maxlength);
    (*system) = choose_option(systemNames, numofsystems, maxlength, "System");
    // Chose System
    if (((*system) > 0) && ((*system) < numofsystems + 1)) {
        // Chose Module
        clear_screen();
        welcome_header(maxlength);
        option_header(toolboxesNames[toolbox - 1], maxlength);
        option_header(systemNames[(*system) - 1], maxlength);
        (*module) = choose_option(moduleNames, numofmodules, maxlength, "Module");
        if (((*module) > 0) && ((*module) < numofmodules + 1)) {
            clear_screen();
            welcome_header(maxlength);
            option_header(toolboxesNames[toolbox - 1], maxlength);
            mixed_header(systemNames[(*system) - 1], moduleNames[(*module) - 1], maxlength);
        }
        else if ((*module) == 0) {
            clear_screen();
            exit(0);
        }
        else {
            clear_screen();
            welcome_header(maxlength);
            option_header(toolboxesNames[toolbox - 1], maxlength);
            option_header(systemNames[(*system) - 1], maxlength);
            invalid_option((*module), "Module", maxlength);
        }
    }
    else if ((*system) == 0) {
        clear_screen();
        exit(0);
    }
    else {
        clear_screen();
        welcome_header(maxlength);
        option_header(toolboxesNames[toolbox - 1], maxlength);
        invalid_option((*system), "System", maxlength);
    }
}