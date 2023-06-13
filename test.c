#include <stdio.h>
#include <stdbool.h>
#include "src/libs/basic.h"


int main(void) {
    char *str = "1\t";

    bool isnum = check_if_string_is_number(str, "int");
    printf("isnum = %s\n", isnum ? "true" : "false");

}