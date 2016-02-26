CC= /usr/bin/gcc
CFLAGS= -c -Iinc -D linux 

OBJS= minimo.o 

.o: .c

libminimo: $(OBJS)
	ar ruv libminimo.a $(OBJS)

clean:
	rm libminimo.a *.o
