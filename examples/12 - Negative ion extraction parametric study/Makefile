CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

simu: simu.o
	$(CC) -o simu simu.o $(LDFLAGS)
simu.o: simu.cpp

clean:
	$(RM) *.o simu


