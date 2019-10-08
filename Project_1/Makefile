    CC     = gcc -g

    CFLAGS = -std=c11
    LFLAGS = -lm

      PROG = pdeSolver
      OBJS = project_1.o \
	  		 utils.o \
             partialDifferential.o

.PHONY: limpa faxina clean purge distclean all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS) $(PROG).o
	$(CC) -o $@ $^ $(LFLAGS)

doc:configura.doxygen $(files)
	doxygen configura.doxygen

limpa:
	@rm -f *~ *.bak

faxina:   clean
	@rm -f *.o core a.out
	@rm -f $(PROG)
