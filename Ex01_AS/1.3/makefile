main.exe: main.cpp ../../RandomGen/random.o
	g++ -Wall -o $@ $^
%.o: %.cpp %.h
	g++ -Wall -c -o $@ $<
clean:
	rm *.o
	rm main.exe
	rm *.out
esegui:
	./main.exe