#
# Makefile for DDFW
#

CC = gcc
CFLAGS = -g -O3 -Wall -Werror -Wno-unused-function -Wno-unused-parameter -std=c99

FILES = ddfw.o clause.o logger.o cnf_parser.o xmalloc.o

ddfw: ddfw.o logger.o clause.o cnf_parser.o xmalloc.o
	$(CC) $(CFLAGS) -o ddfw $(FILES)

verify: Dddfw.o logger.o clause.o cnf_parser.o xmalloc.o
	$(CC) $(CFLAGS) -DDEBUG -o ddfw_verify $(FILES)

ddfw.o: ddfw.c logger.o cnf_parser.o clause.o xmalloc.o
logger.o: logger.c logger.h clause.o
clause.o: clause.c clause.h xmalloc.o
cnf_parser.o: cnf_parser.c cnf_parser.h clause.o
xmalloc.o: xmalloc.c xmalloc.h

Dddfw.o: ddfw.o
	$(CC) $(CFLAGS) -DDEBUG -c -o ddfw.o ddfw.c

clean:
	rm -rf *.o
	rm -f ddfw 
	rm -f ddfw_verify
