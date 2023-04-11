CC = gcc
FLAGS = -Wall -std=c99 -O0 -g

all: 
	$(CC) -o raytracer main.c $(FLAGS)
	raytracer.exe

.PHONY: clean

clean:
	rm -f *.exe *.ppm