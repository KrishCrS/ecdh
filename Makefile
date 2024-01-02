CC = gcc
CFLAGS = -I.
DEPS = ECC.h
OBJ = functions.o ecdh.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

output: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

run: output
	./output

clean:
	rm -f output $(OBJ)
