CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

ParallelGen.exe : ParallelGen.o ../../RandomGen/random.o classi.o
	$(CC) ../../RandomGen/random.o ParallelGen.o classi.o -o ParallelGen.exe
ParallelGen.o : ParallelGen.cpp ParallelGen.h
	$(CC) -c ParallelGen.cpp -o ParallelGen.o $(CFLAGS)
classi.o : classi.cpp classi.h
	$(CC) -c classi.cpp -o classi.o $(CFLAGS)
clean :
	rm *.o ParallelGen.exe
	rm *.out
	rm seed.out
esegui :
	mpiexec -np 4 ./ParallelGen.exe