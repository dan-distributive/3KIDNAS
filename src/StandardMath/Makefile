mFlags=../makeflags  
include $(mFlags)

all: $(StandardMathObj) install

BasicConstants.o:
	$(MakeFObj) BasicConstants.f

random.o:
	$(MakeFObj) random.f

FullCircTrig.o: BasicConstants.o
	$(MakeFObj) FullCircTrig.f

# fdlibm port -- see the note in ../ObjectLists. Compiled with the C
# compiler, not $(MakeFObj); explicit rules (not the *.c pattern rule in
# ../makeflags, which globs all .c files at once) so each stays independent.
fdlibm_sin.o:
	$(CC) $(CFLAGS) -c fdlibm_sin.c

fdlibm_cos.o:
	$(CC) $(CFLAGS) -c fdlibm_cos.c

fdlibm_k_sin.o:
	$(CC) $(CFLAGS) -c fdlibm_k_sin.c

fdlibm_k_cos.o:
	$(CC) $(CFLAGS) -c fdlibm_k_cos.c

fdlibm_rem_pio2_stub.o:
	$(CC) $(CFLAGS) -c fdlibm_rem_pio2_stub.c

fdlibm_log.o:
	$(CC) $(CFLAGS) -c fdlibm_log.c

fdlibm_log1p.o:
	$(CC) $(CFLAGS) -c fdlibm_log1p.c

fdlibm_atanh.o:
	$(CC) $(CFLAGS) -c fdlibm_atanh.c

fdlibm_atan.o:
	$(CC) $(CFLAGS) -c fdlibm_atan.c

fdlibm_exp.o:
	$(CC) $(CFLAGS) -c fdlibm_exp.c

fdlibm_fortran_wrappers.o:
	$(CC) $(CFLAGS) -c fdlibm_fortran_wrappers.c


clean:
	rm -f *.mod;
	rm -f *.o;

install:
	mv *.mod $(LocModLoc);
	mv *.o $(LocBinLoc);
