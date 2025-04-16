#CFLAGS +=  -g  -Wall -mavx512f -mavx512vl -mavx512dq -march=native -fomit-frame-pointer
CFLAGS +=  -O3  -Wall -mavx2 -march=native -fomit-frame-pointer
CC = gcc 

LIB =

LIB_SOURCES = modular_avx256.c cspass_srl3.c cspass_sr.c cspass_rad8.c cspass_rad4.c \
              cs0.c cs1.c cs2.c cs3.c cs4.c cs5.c cs6.c cs7.c cs8.c cs9.c cs10.c cs11.c cs12.c \
              c4_vert.c cs8_vert.c cs16_vert.c ts16_V4.c \
              d0.c d1.c d2.c d3.c d4.c d5.c d6.c d7.c d8.c d9.c d10.c d11.c d12.c \
              uspass_srl3.c uspass_sr.c uspass_rad8.c uspass_rad4.c \
              us0.c us1.c us2.c us3.c us4.c us5.c us6.c us7.c us8.c us9.c us10.c us11.c us12.c\
              u4_vert.c us8_vert.c us16_vert.c uts16_V4.c\
	      rpass_srl3.c rpass_sr.c r0.c r1.c r2.c r3.c r4.c r5.c r6.c r7.c r8.c r9.c r10.c r11.c r12.c\
              vpass_srl3.c vpass_sr.c v0.c v1.c v2.c v3.c v4.c v5.c v6.c v7.c v8.c v9.c v10.c v11.c v12.c mc.c mr.c sc.c sr.c	rpass_srl3_half.c

LIB_OBJECTS =  $(LIB_SOURCES:.c=.o)

UNITTESTS_LIST = $(wildcard test/test_*.c)
UNITTESTS= $(UNITTESTS_LIST:.c=.test)

.SUFFIXES: .c .test

.PHONY: all check clean

all: libntt.a


libntt.a: $(LIB_OBJECTS) 
	$(AR) rc $@ $^
	ranlib $@

check: libntt.a $(UNITTESTS) 
	$(info $(UNITTESTS))
	for t in $(UNITTESTS); do ./$$t ; done 

clean:
	rm libntt.a
	rm -f *.o *~
	rm -f test/*.test test/*.o test/*~ 

.c.o:   $< 
	$(CC) $(CFLAGS) -c $< -o $@

.c.test: $< 
	$(CC) $(CFLAGS) -g -o $@ $< -I. -L. -lntt $(LIB) -lm
