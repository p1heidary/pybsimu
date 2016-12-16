CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

fit: fit.o
	$(CC) -o fit fit.o $(LDFLAGS)
fit.o: fit.cpp

clean:
	$(RM) *.o fit *~

