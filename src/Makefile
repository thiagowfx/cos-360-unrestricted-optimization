CC=g++
CFLAGS=-std=c++11 -g -O2 -Wall
SOURCES=main.cpp lib.hpp
FILES=$(SOURCES) Makefile
EXECUTABLE=main

all: $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -o $(EXECUTABLE)

dist:
	tar zcvf otim-perrotta-$(shell date '+%b-%d-%H-%M').tar.gz $(FILES)
