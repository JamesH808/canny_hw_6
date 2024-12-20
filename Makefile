# Clean Makefile for CENG381
# 10/24/23 EA  adjusted for CENG 381, Fall '23
# Adapted from UCI (RC, ZC, etc.)

SYSTEMC = /opt/pkg/systemc-2.3.4_c++11

INCLUDE = -I. -I$(SYSTEMC)/include
LIBRARY = $(SYSTEMC)/lib-linux64
SYSC_CFLAG = $(INCLUDE) -std=c++11 -L$(LIBRARY) -Xlinker -R -Xlinker $(LIBRARY) -lsystemc -O2

CC = g++
RM = rm -f

VIDEO	= Keck
FRAMES	= $(VIDEO)[0-9][0-9][0-9]_edges.pgm

EXE =	Canny

all: $(EXE)


clean:
	$(RM) $(EXE)
	$(RM) *.bak *~
	$(RM) *.o *.ti gmon.out
	$(RM) $(IMG_OUT)
	$(RM) $(FRAMES)

cleanall: clean
	$(RM) *.log

	
# Assignment 6 
# Test bench model in SystemC

Canny: Canny.cpp
	$(CC) $(SYSC_CFLAG) -pg $< -o $@

test: Canny video/$(FRAMES)
	ulimit -s 128000; ./$<
	set -e; \
	for f in video/$(FRAMES); do \
		diff `basename $$f` $$f; \
	done

# EOF
