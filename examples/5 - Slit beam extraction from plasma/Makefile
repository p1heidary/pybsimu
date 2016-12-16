CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

slit: slit.o
	$(CC) -o slit slit.o $(LDFLAGS)
slit.o: slit.cpp

clean:
	$(RM) *.o slit 

