#
# Makefile for DDFW
#

CC = gcc
CFLAGS = -g -O3 -Wall -Werror -Wno-unused-function -Wno-unused-parameter -std=c99

FILES = main.o ddfw.o logger.o clause.o neighborhood.o weight_transfer.o cnf_parser.o verifier.o xmalloc.o

ddfw: main.o ddfw.o logger.o clause.o neighborhood.o weight_transfer.o cnf_parser.o verifier.o xmalloc.o
	$(CC) $(CFLAGS) -o ddfw $(FILES)

verify: main.o Dddfw.o logger.o clause.o neighborhood.o weight_transfer.o cnf_parser.o Dverifier.o xmalloc.o
	$(CC) $(CFLAGS) -DDEBUG -o ddfw_verify $(FILES)

main.o: main.c
ddfw.o: ddfw.c ddfw.h
logger.o: logger.c logger.h
clause.o: clause.c clause.h
neighborhood.o: neighborhood.c neighborhood.h
weight_transfer.o: weight_transfer.c weight_transfer.h
cnf_parser.o: cnf_parser.c cnf_parser.h
verifier.o: verifier.c verifier.h
xmalloc.o: xmalloc.c xmalloc.h

Dddfw.o: ddfw.o
	$(CC) $(CFLAGS) -DDEBUG -c -o ddfw.o ddfw.c
Dverifier.o: verifier.o
	$(CC) $(CFLAGS) -DDEBUG -c -o verifier.o verifier.c

clean:
	rm -rf *.o
	rm -f ddfw 
	rm -f ddfw_verify
