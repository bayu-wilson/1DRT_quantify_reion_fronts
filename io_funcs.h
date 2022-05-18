#ifndef IO_FUNCS_H
#define IO_FUNCS_H

void get_j_lya();
void read_grid();
void read_source();
void make_output();
void write_otf();
void write_otf_bayu();
void write_spectrum();
void write_gas(bool bool_initial_gas);
void loading_bar(double number, double max_number);

#endif
