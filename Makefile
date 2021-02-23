#
# Makefile for DDFW
#

CC = gcc
CFLAGS = -g -O3 -Wall -Werror -Wno-unused-function -Wno-unused-parameter -std=c99

FILES = main.o ddfw.o logger.o clause.o neighborhood.o weight_transfer.o cnf_parser.o xmalloc.o

ddfw: main.o ddfw.o logger.o clause.o neighborhood.o weight_transfer.o cnf_parser.o xmalloc.o
	$(CC) $(CFLAGS) -o ddfw $(FILES)

verify: main.o Dddfw.o logger.o clause.o neighborhood.o weight_transfer.o cnf_parser.o xmalloc.o
	$(CC) $(CFLAGS) -DDEBUG -o ddfw_verify $(FILES)

main.o: main.c ddfw.o logger.o cnf_parser.o clause.o xmalloc.o
ddfw.o: ddfw.c logger.o cnf_parser.o clause.o xmalloc.o
logger.o: logger.c logger.h clause.o
clause.o: clause.c clause.h xmalloc.o
neighborhood.o: neighborhood.c neighborhood.h
weight_transfer.o: weight_transfer.c weight_transfer.h
cnf_parser.o: cnf_parser.c cnf_parser.h clause.o
xmalloc.o: xmalloc.c xmalloc.h

Dddfw.o: ddfw.o
	$(CC) $(CFLAGS) -DDEBUG -c -o ddfw.o ddfw.c

clean:
	rm -rf *.o
	rm -f ddfw 
	rm -f ddfw_verify
