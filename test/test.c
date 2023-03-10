#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

void partition(int mode, size_t maxlength) {
    if (mode == 1) {
        printf(" ");
        for (size_t i = 0; i < maxlength-2; i++) {
            printf("=");
        }
        printf(" ");
    }
    else if (mode == 2) {
        printf(" ");
        for (size_t i = 0; i < maxlength-2; i++) {
            printf("-");
        }
        printf(" ");
    }
    else {
        printf("DEBUG WARNING: Please select 1 for main partition or 2 for secondary partition for the interface programming");
        exit(0);
    }
    printf("\n");
}

int int_length(int value) {
    if (value == 0) {
        return 1;
    } 
    else {
        return (floor(log10(abs(value))) + 1);
    }
}

void write_sys_parameters(int npar, double *par, size_t maxlength, double percname) {
    int spcname = maxlength - (1 - percname)*maxlength; // Space for name of printer variable
    int spcvalue = maxlength - percname*maxlength;      // Space for value of variable
    partition(2, maxlength);
    printf("  System Parameters\n");
    partition(2, maxlength);
    for (int i = 0; i < npar; i++) {
        //printf("%s%d%-*s %-*g\n", "  par[", i, spcname - 7, "]:", spcvalue, par[i]);
        //printf("%s%d%-*s %-*g\n", "  par[", i, spcname - 7, "]:", spcvalue, par[i]);
        printf("%s%d%-*s %-*g\n", "  par[", i, spcname - 7 - int_length(i), "]:", spcvalue, par[i]);
    }
}

int store_indexes_strlen(double n, int *indexes) {
    int slen = 0;
    // determine the length necessary to store the string containing the list of values of the indexes with comma and space(, ) 
    for (int i = 0; i < n; i++) {
        if (i == n - 1) {
            // previous slen + new int length + length of the '\0' char
            slen = slen + int_length(indexes[i]) + 1;
        } else {
            // previous slen + new int length + length of the comma + length of the space
            slen = slen + int_length(indexes[i]) + 2;
        }        
    }

    return slen;
}

char* convert_list_of_int_to_string_working(int slen, int n, int *indexes, int spcvalue) {
    // Allocate memory to handle operations
    char *string = malloc(slen * sizeof(char));  // Final String
    char *buffer = malloc(slen * sizeof(char));  // Buffer to store each index at a time
    char *newrow = malloc(spcvalue * sizeof(char)); // String to store each row of the final string
    // Sweep across all indexes
    for (int i = 0; i < n; i++) {
        // Check if is the last index
        if (i == n - 1) {
            // Put the final index into the buffer
            sprintf(buffer, "%d", indexes[i]);
            // Check if the length of buffer and newrow together are bigger than the maximum space for values
            if (strlen(newrow) + strlen(buffer) > spcvalue) {
                // Save in the final string the new row
                strcat(string, newrow);
                // Empty the row string
                strcpy(newrow, "");
                // Add the final index in the buffer in the new line
                sprintf(buffer, "\n%d", indexes[i]);
            }
            // Update newrow with the new value
            strcat(newrow, buffer);
            // Update final string with the last row
            strcat(string, newrow);
        }
        else {
            // Put the next index into the buffer
            sprintf(buffer, "%d, ", indexes[i]);
            // Check if the length of buffer and newrow together are bigger than the maximum space for values
            if (strlen(newrow) + strlen(buffer) > spcvalue) {
                // Save in the final string the new row
                strcat(string, newrow);
                // Empty the row string
                strcpy(newrow, "");
                // Add the next index in the buffer into the new line
                sprintf(buffer, "\n%d, ", indexes[i]);
            }
            // Update newrow with the new value
            strcat(newrow, buffer);
        }
    }
    // Free Memory
    free(newrow); free(buffer);
    // Return the final string (the programmer is responsible to free the final string after the function call)
    return string;
}

void print_list_of_indexes(int n, int *indexes, int spcvalue, int spcname, char* name) {
    // Allocate memory to handle operations
    char *buffer = malloc((spcvalue + 1) * sizeof(char));  // Buffer to store each index at a time
    char *newrow = malloc((spcvalue + 1) * sizeof(char)); // String to store each row
    char *formatedrow = malloc((spcvalue + spcname + 2) * sizeof(char)); // String to store each formatted row
    // Declare counter for number of rows found
    int rows = 0; 
    // Sweep across all indexes
    for (int i = 0; i < n; i++) {
        // Check if is the last index
        if (i == n - 1) {
            // Put the final index into the buffer
            sprintf(buffer, "%d", indexes[i]);
            // Check if the length of buffer and newrow together are bigger than the maximum space for values
            if (strlen(newrow) + strlen(buffer) > spcvalue) {
                // Format the row with names
                sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
                printf("%s\n", formatedrow);
                // Empty the row string
                strcpy(newrow, "");
                strcpy(formatedrow, "");
            }
            // Update newrow with the new value
            strcat(newrow, buffer);
            // Format the row with names
            sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
            printf("%s\n", formatedrow);
        }
        else {
            // Put the next index into the buffer
            sprintf(buffer, "%d, ", indexes[i]);
            // Check if the length of buffer and newrow together are bigger than the maximum space for values
            if (strlen(newrow) + strlen(buffer) > spcvalue) {
                rows = rows + 1;
                // Format the row with names
                if (rows == 1) {
                    sprintf(formatedrow, "%-*s %-*s", spcname, name, spcvalue, newrow);
                }
                else {
                    sprintf(formatedrow, "%-*s %-*s", spcname, "", spcvalue, newrow);
                }
                // Save in the final string the new row
                printf("%s\n", formatedrow);
                // Empty the row string
                strcpy(newrow, "");
                strcpy(formatedrow, "");
            }
            // Update newrow with the new value
            strcat(newrow, buffer);
        }
    }
    // Free Memory
    free(newrow); free(buffer); free(formatedrow);
}

int main(void) {
    // Parameters related to printing information
    size_t maxLen = 71;             // Max length of the info printed on the screen and on info file
    double percName = 0.6;          // Percentage of space occuped by the name of the quantity printed
    int spcvalue = maxLen - percName*maxLen;
    int spcname = maxLen - (1 - percName)*maxLen;  // Space for name of printer variable
    
    //write_sys_parameters(npar, par, maxLen, percName);
    int indexes[] = {0, 5, 10, 15, 20, 21, 22, 23, 105, 110, 114, 256, 789, 800, 998, 760, 150000000};
    int slen = store_indexes_strlen(17, indexes);
    
    print_list_of_indexes(17, indexes, spcvalue, spcname, "  Printed on File Indexes:");

}