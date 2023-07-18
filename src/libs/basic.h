#ifndef BASIC_H
#define BASIC_H

int count_int_digits(int num);
void **alloc_2D_array(int rows, int columns, size_t elementSize);
char **alloc_string_array(int nstrings, size_t maxlen);
void *copy_pointer(const void* src, int size, size_t element_size);
void **copy_2D_pointer(const void** src, int rows, int cols, size_t element_size);
void free_mem(void* first, ...);
void free_2D_mem(void **array_2D, int nrows);
void close_files(int num_files, ...);
void file_safety_check(FILE *file);
void ptr_safety_check(void *ptr, char *ptr_name);
void double_ptr_safety_check(void** ptr, char *ptr_name);
int get_largest_element_int_array(int *arr, int size);

#endif