CC=g++
CFLAGS = -Iinclude -Isrc -static -lstdc++ -Llib/lib-mingw-w64 -O3 -s
LIBS = -lglfw3 -lgdi32 -lopengl32

build/SaccEngine.exe: src/SaccEngine.cpp
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

clean:
	rm -f *.o *~