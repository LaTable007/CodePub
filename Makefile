CC = gcc-13
# librairies de PRIMME
LIBP = -L./primme/ -lprimme -L./EVSL_1.1.1/EVSL/lib/ -levsl
# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/ -I./EVSL_1.1.1/EVSL/INC 
# toutes les librairies
LIB = $(LIBP) -lm -lblas -llapack

COPT = -g -Wall

default: main

clean: 
	rm *.o 
	rm main

main: main.c prob.o time.o interface_primme.o matvec.o interface_evsl.o residu.o prob.h time.h interface_primme.h matvec.h interface_evsl.h residu.h
	$(CC) $(COPT) $^ -o $@ $(LIB)

%.o: %.c
	$(CC) $(COPT) -c $< -o $@ $(INCP)


