CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

einzel3d: einzel3d.o
	$(CC) -o einzel3d einzel3d.o $(LDFLAGS)
einzel3d.o: einzel3d.cpp

clean:
	$(RM) *.o einzel3d 

