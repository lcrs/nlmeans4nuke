CPPFLAGS = -I/usr/local/nuke/include
CFLAGS = -Wall -DUSE_GLEW -c -O3 -march=pentium-m -mfpmath=sse
LFLAGS = -shared 

Nlmeans.so: nlmeans.o nlmeans_lib.o Makefile
	g++ $(LFLAGS) -o Nlmeans.so nlmeans.o nlmeans_lib.o

nlmeans.o: nlmeans.cpp Makefile
	g++ $(CPPFLAGS) $(CFLAGS) -o nlmeans.o nlmeans.cpp

nlmeans_lib.o: nlmeans_lib.cpp nlmeans_lib.h Makefile
	g++ $(CPPFLAGS) $(CFLAGS) -o nlmeans_lib.o nlmeans_lib.cpp

clean:
	rm -f nlmeans.o nlmeans_lib.o Nlmeans.so

install: Nlmeans.so Makefile
	cp Nlmeans.so ~/.nuke/
