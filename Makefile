CC=g++

RM=/bin/rm -f
# INTEL= -mcmodel=medium # temp: from chris
# CFLAGS=$(INTEL) -O3 -std=c++17 # temp: from chris # note: '-fopenmp' needed for parallelization
# Plog is a C++ logging library. see https://github.com/SergiusTheBest/plog for more info
# -I option to specify an alternate include directory
# CFLAGS=-std=c++17 -I/Volumes/Extreme\ SSD/utility/plog-master/include/ -Wall -Xpreprocessor -fopenmp #-O3
CFLAGS=-std=c++17 -Wall -fopenmp -O3
###INC = -I/opt/homebrew/opt/libomp/include/  -I/Volumes/Extreme\ SSD/utility/plog-master/include/ #-I/usr/local/include
###LIB = -L/opt/homebrew/opt/libomp/lib/ -lomp #-L/usr/local/lib -lomp
EXE = 1d_radtransfer.x
OBJ = main.o init_funcs.o general_funcs.o global_variables.o io_funcs.o rates.o rt_funcs.o gas_funcs.o cosmo_funcs.o data_funcs.o
###LOG = Logfile.log

main: $(OBJ)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ)

.cc.o: MAKEFILE
	$(CC) $(CFLAGS) -o $*.o -c $*.cc

clean:
	$(RM) $(EXE) $(OBJ)
