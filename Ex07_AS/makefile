Monte_Carlo_NVT.exe: Monte_Carlo_NVT.cpp ../RandomGen/random.o
	g++ -Wall -o $@ $^
%.o: %.cpp %.h
	g++ -Wall -c -o $@ $<
clean:
	rm *.o
	rm main.exe
	rm seed.out
esegui:
	./Monte_Carlo_NVT.exe
