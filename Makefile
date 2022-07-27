OBJS = options.o word.o seq_overlap.o parse_fasta.o split.o
MOTIF_OBJS = word.o parse_fasta.o mpi_util.o motif_options.o parse_codon.o split.o
	
CC = mpic++
cc = mpicc

# Use -DATAx2 for MMX-based 32 bit alignments of 2 sequences in parallel (-mmmx required)
# Use -DATAx4 for SSE2-based 32 bit alignments of 4 sequences in parallel (-msse4.1 required)
# Use -DATAx8 for AVX2-based 32 bit alignments of 8 sequences in parallel (-mavx2 required)
#
# Use -DEBUG to enable bounds checking
FLAGS = -O3 -Wall -mavx2 -DATA32x4

# Add '-g' to include debug symbols for use with 'perf'
PROFILE = -g#-pg
OPENMP = -fopenmp

INC = -I. -I$(HOME)/gsl/include
LIBS = -lm -lz $(HOME)/gsl/lib/libgsl.a

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(PROFILE) $(OPENMP) $(INC) -c $<
 
.c.o:
	$(cc) $(FLAGS) $(PROFILE) $(OPENMP) $(INC) -c $<

all: genome-diff motif

motif: $(MOTIF_OBJS) motif.o
	$(CC) $(PROFILE) -o motif $(MOTIF_OBJS) $(OPENMP) motif.o $(LIBS)

arborate: $(OBJS) arborate.o
	$(CC) $(PROFILE) -o arborate $(OBJS) $(OPENMP) arborate.o $(LIBS)

genome-diff: $(OBJS) genome_diff.o
	$(CC) $(PROFILE) -o genome-diff $(OBJS) $(OPENMP) genome_diff.o $(LIBS)

clean:
	-rm -f *.o


