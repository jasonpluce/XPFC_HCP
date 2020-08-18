
OBJS = main.o GB_system.o get_size.o make_C.o plan_fftw.o set_IC.o evolve.o test_strain.o write_out.o visit_writer.o

CC = g++ -O3
C = gcc -c
DEBUG = -g
CFLAGS = -Wall -c -I${FFTW_INCLUDE} -I/home/jpluce/blitz-0.10 -std=c++11
#CFLAGS = -Wall -c -I${FFTW_INCLUDE} -I/home/jpluce/blitz-0.10 -std=c++11 -fopenmp

LFLAGS = -Wall
LIBS = -L${FFTW_LIB} -fopenmp -lfftw3 -lfftw3_omp -lm

PFC : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o PFC $(LIBS)
	
main.o: header.h main.cpp
	$(CC) $(CFLAGS) main.cpp
	
GB_system.o: header.h GB_system.cpp
	$(CC) $(CFLAGS) GB_system.cpp
	
get_size.o: header.h get_size.cpp
	$(CC) $(CFLAGS) get_size.cpp

make_C.o: header.h make_C.cpp
	$(CC) $(CFLAGS) make_C.cpp

plan_fftw.o: header.h plan_fftw.cpp
	$(CC) $(CFLAGS) plan_fftw.cpp

set_IC.o: header.h set_IC.cpp
	$(CC) $(CFLAGS) set_IC.cpp

evolve.o: header.h evolve.cpp
	$(CC) $(CFLAGS) evolve.cpp

write_out.o: header.h visit_writer.h write_out.cpp
	$(CC) $(CFLAGS) write_out.cpp

test_strain.o: header.h test_strain.cpp
	$(CC) $(CFLAGS) test_strain.cpp

visit_writer.o: visit_writer.h visit_writer.c
	$(C) visit_writer.c 
	
clean:
	\rm *.o PFC

run: PFC
	#> blank.vtk
	#> blank.txt
	#> PFC.e
	#> PFC.o
	#rm *.vtk
	#rm *.txt
	#rm PFC.e*
	#rm PFC.o*
	./PFC
