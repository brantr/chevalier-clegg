
EXEC    = chevalier

CC      = clang++              # sets the C-compiler

OPTIONS = $(OPTIMIZE) $(OPT)


OBJS   = main.o chevalier.o

INCL   = chevalier.hpp

CFLAGS = $(OPTIONS)

LIBS   =  -lm -lgsl -lgslcblas

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)
