SimulAnn: SimulAnn.cpp ../../RandomGen/random.o classi.o
	g++ -Wall -o $@ $^
%.o: %.cpp %.h
	g++ -Wall --std=c++11 -c -o $@ $<
clean:
	rm *.o
	rm SimulAnn
	rm seed.out
esegui:
	./SimulAnn