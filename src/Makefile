#
# Makefile for DDFW
#

CC = gcc
CFLAGS = -g -O3 -Wall -Werror -Wno-unused-function -Wno-unused-parameter -std=c99

FILES = main.o ddfw.o logger.o assignment.o neighborhood.o weight_reducer.o weight_transfer.o cnf_parser.o verifier.o initializer.o global_data.o xmalloc.o

ddfw: main.o ddfw.o logger.o assignment.o neighborhood.o weight_reducer.o weight_transfer.o cnf_parser.o verifier.o initializer.o global_data.o xmalloc.o
	$(CC) $(CFLAGS) -o ddfw $(FILES)

verify: main.o Dddfw.o logger.o assignment.o neighborhood.o weight_reducer.o weight_transfer.o cnf_parser.o Dverifier.o initializer.o global_data.o xmalloc.o
	$(CC) $(CFLAGS) -DDEBUG -o ddfw_verify $(FILES)

main.o: main.c
ddfw.o: ddfw.c ddfw.h
logger.o: logger.c logger.h
assignment.o: assignment.c assignment.h
neighborhood.o: neighborhood.c neighborhood.h
weight_reducer.o: weight_reducer.c weight_reducer.h
weight_transfer.o: weight_transfer.c weight_transfer.h
cnf_parser.o: cnf_parser.c cnf_parser.h
verifier.o: verifier.c verifier.h
initializer.o: initializer.c initializer.h
global_data.o: global_data.c global_data.h
xmalloc.o: xmalloc.c xmalloc.h

Dddfw.o: ddfw.o
	$(CC) $(CFLAGS) -DDEBUG -c -o ddfw.o ddfw.c
Dverifier.o: verifier.o
	$(CC) $(CFLAGS) -DDEBUG -c -o verifier.o verifier.c

clean:
	rm -rf *.o
	rm -f ddfw 
	rm -f ddfw_verify
