#ifndef BASIC_H
#define BASIC_H

int count_int_digits(int num);
void free_mem(void* first, ...);
void free_2D_mem(void **array_2D, int nrows);
void close_files(int num_files, ...);
void file_safety_check(FILE *file);
void ptr_safety_check(void *ptr, char *ptr_name);

#endif