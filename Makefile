CC = gcc-13
# librairies de PRIMME
LIBP = -L./primme/ -lprimme
# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/ 
# toutes les librairies
LIB = $(LIBP) -lm -lblas -llapack

COPT = -g -Wall

default: main

clean: 
	rm *.o 
	rm main

main: main.c prob.o time.o interface_primme.o
	$(CC) $(COPT) $^ -o $@ $(LIB)

%.o: %.c
	$(CC) $(COPT) -c $< -o $@ $(INCP)


