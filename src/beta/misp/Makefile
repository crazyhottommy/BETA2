SOURCES = motif.c misp.c
CC = cc
LIB = -lz -lm


all: kseq.h misp.c motif.c
	$(CC) -Wall -g $(SOURCES) -o misp -O3 $(LIB)

dbg: kseq.h misp.c motif.c
	$(CC) -g $(SOURCES) -o misp $(LIB)

clean:
	rm -f *.o misp
	rm test_*
	rm test2_*

check-syntax:
	$(CC) -Wall -Wextra -pedantic -fsyntax-only $(SOURCES) $(LIB)


test:
	make all
	time ./misp test.seq database/cistrome.db 0.001 all test
	time ./misp test.seq database/cistrome.db 0.001 EN0055 test
	time ./misp test2.seq database/cistrome.db 0.001 all test2
	time ./misp test2.seq database/cistrome.db 0.001 all test2
	time ./misp test3.seq database/cistrome.db 0.001 hPDI111 test3
	time ./misp test3.seq database/cistrome.db 0.001 M01195 test3
