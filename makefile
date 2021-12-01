all: aout

aout: mod_spectra_heom.o spectra_heom.o
	ifort -o aout -O3 -qopenmp mod_spectra_heom.o spectra_heom.o

%.o: %.f90
	ifort -c -qopenmp $<

clean:
	rm *.o aout

