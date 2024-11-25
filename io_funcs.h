#ifndef IO_FUNCS_H
#define IO_FUNCS_H

void get_j_lya();
// void read_grid();
void read_grid_mmc();
void read_source();
void make_output();
void write_otf();
void write_otf_bayu();
void write_otf_fred(char output_string[]);
void write_otf_spectrum(char output_string[]);
void write_spectrum();
void write_gas(char output_path[]);
void loading_bar(double number, double max_number);
void write_Inu(char output_string[]); 

#endif
