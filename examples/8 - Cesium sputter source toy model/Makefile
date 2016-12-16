CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

sputter: sputter.o
	$(CC) -o sputter sputter.o $(LDFLAGS)
sputter.o: sputter.cpp

clean:
	$(RM) *.o sputter 

