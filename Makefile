CC     = gcc
CFLAGS = -Wall -O2
LIBS   = -lz

TARGET = bntools
MODULES = $(patsubst src/%.c,%,$(wildcard src/*.c))

.PHONY: all clean

all: ${TARGET}

clean:
	@rm -vrf tmp/ ${TARGET}

${TARGET}: ${MODULES:%=tmp/%.o}
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

tmp/%.o: src/%.c
	@[ -d ${@D} ] || mkdir -pv ${@D}
	${CC} -c ${CFLAGS} -o $@ $^
