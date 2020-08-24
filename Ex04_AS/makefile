MolDyn_NVE.exe: MolDyn_NVE.cpp
	g++ -Wall -o $@ $^
%.o: %.cpp %.h
	g++ -Wall -c -o $@ $<
esegui:
	./MolDyn_NVE.exe