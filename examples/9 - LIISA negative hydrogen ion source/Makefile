CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

liisa: liisa.o
	$(CC) -o liisa liisa.o $(LDFLAGS)
liisa.o: liisa.cpp

liisa_cont: liisa_cont.o
	$(CC) -o liisa_cont liisa_cont.o $(LDFLAGS)
liisa_cont.o: liisa_cont.cpp

liisa_mg: liisa_mg.o
	$(CC) -o liisa_mg liisa_mg.o $(LDFLAGS)
liisa_mg.o: liisa_mg.cpp

clean:
	$(RM) *.o liisa liisa_cont liisa_mg *~

