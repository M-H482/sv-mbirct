# Change the compiler (e.g. icc,gcc,clang) in the following variable,
# or provide an argument as "make CC=gcc". If provided, the argument will
# cause the following CC variable definition to be ignored.
CC = gcc

# common flags, libraries
CFLAGS = -std=c11 -O3
CFLAGS += -Wall
LIBS = -lm

ifeq ($(CC),icc)
  CFLAGS += -qopenmp
  CFLAGS += -no-prec-div -restrict -ipo -inline-calloc -qopt-calloc -no-ansi-alias
  CFLAGS += -xCORE-AVX2
  #CFLAGS += -qopenmp-link=static
  #CFLAGS += -static-intel -wd10237
  #CFLAGS += -prof-gen
  #CFLAGS += -prof-use
  #CFLAGS += -no-vec -no-simd
  #CFLAGS += -qopt-report=4 -qopt-report-phase=vec
  #CFLAGS += -w3
  #CFLAGS += -debug

else ifeq ($(CC),clang)
  CFLAGS += -Xclang -fopenmp
  CFLAGS += -Wno-unknown-pragmas
  LIBS += -lomp

else ifeq ($(CC),gcc)
  CFLAGS += -fopenmp
  CFLAGS += -Wno-unknown-pragmas
else
  CFLAGS += -fopenmp
endif


all: mbir_ct mbir_ct_stream clean

clean:
	rm *.o

# Common Objects (Full Reconstruction)
OBJ = initialize.o recon3d.o heap.o icd3d.o allocate.o MBIRModularUtils.o

# Stream Objects (Matrix Generation Only - NO recon3d/icd3d)
OBJ_STREAM = initialize.o allocate.o MBIRModularUtils.o

# Original Version (Complete Reconstruction)
mbir_ct: mbir_ct.o A_comp.o $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)
	mv $@ ../bin

# Streaming Version (Matrix Generation Only)
mbir_ct_stream: mbir_ct_stream.o A_comp_stream.o $(OBJ_STREAM)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)
	mv $@ ../bin

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

