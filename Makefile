# Copyright (C) 2013 Lile Wang
# For full notice, see "main.cpp" and "COPYING".

CC = g++
Option = -O3 -fopenmp -Wall

Objs = main.o kdtree.o correlate.o read_data.o input.o \
	driver.o parallel.o

Out_file = kdtpcf

main : $(Objs)
	$(CC) $(Option) $(Objs) -o $(Out_file)

main.o : main.cpp
	$(CC) $(Option) -c -o main.o main.cpp

kdtree.o : kdtree.cpp kdtree.h
	$(CC) $(Option) -c -o kdtree.o kdtree.cpp

correlate.o : correlate.cpp correlate.h
	$(CC) $(Option) -c -o correlate.o correlate.cpp

read_data.o : read_data.cpp read_data.h
	$(CC) $(Option) -c -o read_data.o read_data.cpp

input.o : input.cpp input.h
	$(CC) $(Option) -c -o input.o input.cpp

driver.o : driver.cpp driver.h
	$(CC) $(Option) -c -o driver.o driver.cpp

parallel.o : parallel.cpp parallel.h
	$(CC) $(Option) -c -o parallel.o parallel.cpp

clean:
	rm *.o *~ kdtpcf
