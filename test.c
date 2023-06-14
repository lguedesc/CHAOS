#include <stdio.h>
#include <stdbool.h>
#include "src/libs/basic.h"

bool isNegativeNumber(const char* str) {
    if (str[0] == '-') {  // Check if the first character is a minus sign
        // Check if the remaining characters are digits
        for (int i = 1; str[i] != '\0'; i++) {
            if (str[i] < '0' || str[i] > '9') {
                return false;  // Found a non-digit character
            }
        }
        return true;  // All characters are digits after the minus sign
    }
    return false;  // First character is not a minus sign
}

int main(void) {
    char *str = "-1\n\t";

    bool isnum = check_if_string_is_number(str, "double");
    bool isneg = isNegativeNumber(str);
    printf("isnum = %s\n", isnum ? "true" : "false");
    printf("isneg = %s\n", isnum ? "true" : "false");

}