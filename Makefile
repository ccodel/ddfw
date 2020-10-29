#
# Makefile for DDFW
#

CC = gcc
CFLAGS = -g -O2 -Wall -Werror -Wno-unused-function -Wno-unused-parameter

FILES = ddfw.o clause.o logger.o cnf_parser.o

ddfw: ddfw.o logger.o clause.o cnf_parser.o
	$(CC) $(CFLAGS) -o ddfw $(FILES)

ddfw.o: ddfw.c logger.o cnf_parser.o clause.o 
logger.o: logger.c logger.h clause.o
clause.o: clause.c clause.h 
cnf_parser.o: cnf_parser.c cnf_parser.h clause.o

clean:
	rm -rf *.o
	rm -f ddfw 
