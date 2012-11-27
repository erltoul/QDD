COMPILERFLAGS = -O3
LIBS = -lfftw3 -lm


EXEC = create_wisdom

OBJS = create_wisdom.o fftw.o

%.o: %.mod

$(EXEC): $(OBJS)
	ifort $(COMPILERFLAGS) -o $@ $(OBJS) $(LIBS)

create_wisdom.o: create_wisdom.F90 define.h makefile fftw.o
	ifort $(COMPILERFLAGS) -c $<

fftw.o: fftw.F90 makefile
	ifort $(COMPILERFLAGS) -c $<

clean:
	rm -rf *.o *.mod
