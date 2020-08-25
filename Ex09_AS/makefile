Genetic.exe: Genetic.cpp ../RandomGen/random.o classi.o
	g++ -Wall -o $@ $^
%.o: %.cpp %.h
	g++ -Wall --std=c++11 -c -o $@ $<
clean:
	rm *.o
	rm main.exe
	rm seed.out
esegui:
	./Genetic.exe