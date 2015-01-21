CC     = gcc
CFLAGS = -Wall
LIBS   = -lz

ifeq ("${DEBUG}", "")
CFLAGS += -O2
else
CFLAGS += -g
endif

APP_VER := $(shell src/version.sh)

TARGET = bntools
MODULES = $(patsubst src/%.c,%,$(wildcard src/*.c))

.PHONY: all clean

all: ${TARGET}

clean:
	@rm -vrf tmp/ ${TARGET} src/version.h

${TARGET}: ${MODULES:%=tmp/%.o}
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

tmp/%.o: src/%.c
	${CC} -c ${CFLAGS} -o $@ $<

sinclude ${MODULES:%=tmp/%.d}
tmp/%.d: src/%.c
	@echo "Parsing dependency for '$<'"
	@[ -d ${@D} ] || mkdir -pv ${@D}
	@${CC} -MM $< -MT ${@:%.d=%.o} | sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $@
