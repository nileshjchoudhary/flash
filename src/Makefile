
FC = fc
all: $(FC)

FC_SRC = $(wildcard *.c)
FC_DIR = $(notdir $(FC_SRC))
FC_OBJ = $(patsubst %.c,%.o,$(FC_DIR))

LIBS = -lm
CC = mpicc

$(FC): $(FC_OBJ) 
	${CC} -o $(FC) $(FC_OBJ) $(LDFLAGS) $(LIBS) 

clean:
	@$(RM) -fv *.o $(FC)
