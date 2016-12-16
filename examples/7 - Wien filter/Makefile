CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

wien: wien.o
	$(CC) -o wien wien.o $(LDFLAGS)
wien.o: wien.cpp

clean:
	$(RM) *.o wien 

