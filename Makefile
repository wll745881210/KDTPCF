# Copyright (C) 2013 Lile Wang
# For full notice, see "main.cpp" and "COPYING".

CC = g++
Option = -O3 -Wall

Objs = main.o kdtree.o correlate.o read_data.o input.o driver.o

Out_file = test_corr

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

