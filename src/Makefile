
FC = fc
all: $(FC)



FC_SRC = $(wildcard *.c)
FC_DIR = $(notdir $(FC_SRC))
FC_OBJ = $(patsubst %.c,%.o,$(FC_DIR))

LIBS = -lm
CC = mpicc
CFLAGS = -g -O2 -Wall -Wextra  -Wmissing-declarations -Wmissing-prototypes -Wimplicit-function-declaration -Wunused -Wuninitialized 

$(FC): $(FC_OBJ) 
	${CC} -O2 -g -o $(FC) $(FC_OBJ) $(LDFLAGS) $(LIBS) 

clean:
	@$(RM) -fv *.o $(FC)
