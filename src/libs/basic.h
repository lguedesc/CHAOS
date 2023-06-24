#include <stdbool.h>

int count_int_digits(int num);
void free_mem(void* first, ...);
void free_2D_mem(void **array_2D, int nrows);
bool check_if_string_is_number_old(const char* str);
bool check_if_string_is_number(const char *str, const char *type, bool only_positive);
bool check_if_string_is_negative_number(const char *str);
char **malloc_string_array(int nstrings, int maxlen);