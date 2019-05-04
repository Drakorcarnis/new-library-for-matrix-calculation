#
# Coverage and profiling data
#
COV_DIR = ../cov
PROF_DIR = ../prof
#
# Common compiler flags
#
CC     = gcc -pipe
CFLAGS = -Wall -Werror -Wextra
LDFLAGS = -lm -lgomp
SOFLAGS = -shared -fPIC -fopenmp -march=native

BINDIR = ../bin
BINDBGDIR = $(BINDIR)/debug
BINRELDIR = $(BINDIR)/release

OBJDIR = ../obj
OBJDBGDIR = $(OBJDIR)/debug
OBJRELDIR = $(OBJDIR)/release

LIBDIR = ../lib
LIBDBGDIR = $(LIBDIR)/debug
LIBRELDIR = $(LIBDIR)/release

OBJLIBDIR = $(OBJDIR)/lib
OBJLIBDBGDIR = $(OBJLIBDIR)/debug
OBJLIBRELDIR = $(OBJLIBDIR)/release
#
# Project files
#
INCLUDES = includes
LIB_SRCS = matrix.c matrix_tools.c plu.c cholesky.c check.c
TEST_SRCS = test.c
REG_SRCS = regression.c
LIB_OBJS = $(LIB_SRCS:.c=.o)
TEST_OBJS = $(TEST_SRCS:.c=.o)
REG_OBJS = $(REG_SRCS:.c=.o)
TEST_EXE  = test
REGRESSION_EXE = regression
LIB  = matrix

#
# Lib Debug build settings
#
LIBDBGEXE = $(LIBDBGDIR)/lib$(LIB).so
LIBDBGOBJS = $(addprefix $(OBJLIBDBGDIR)/, $(LIB_OBJS))
DBGCFLAGS = -g -O0 -DDEBUG --coverage
DBGLDFLAGS = -L$(LIBDBGDIR) -l$(LIB)
$(OBJLIBDBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(SOFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES)
$(OBJDBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES)
#
# Lib Release build settings
#
LIBRELEXE = $(LIBRELDIR)/lib$(LIB).so
LIBRELOBJS = $(addprefix $(OBJLIBRELDIR)/, $(LIB_OBJS))
RELCFLAGS = -Ofast -DNDEBUG  -fno-trapping-math -fopt-info-vec-missed 
RELLDFLAGS = -L$(LIBRELDIR) -l$(LIB)
$(OBJLIBRELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(SOFLAGS) $(RELCFLAGS) -o $@ $< -I$(INCLUDES)
$(OBJRELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(RELCFLAGS) -o $@ $< -I$(INCLUDES)
#
# Test Debug build settings
#
TESTDBGEXE = $(BINDBGDIR)/$(TEST_EXE)
TESTDBGOBJS = $(addprefix $(OBJDBGDIR)/, $(TEST_OBJS))

#
# Test Release build settings
#
TESTRELEXE = $(BINRELDIR)/$(TEST_EXE)
TESTRELOBJS = $(addprefix $(OBJRELDIR)/, $(TEST_OBJS))

#
# Regression build settings
#
REGDBGEXE = $(BINDBGDIR)/$(REGRESSION_EXE)
REGDBGOBJS = $(addprefix $(OBJDBGDIR)/, $(REG_OBJS))

REGRELEXE = $(BINRELDIR)/$(REGRESSION_EXE)
REGRELOBJS = $(addprefix $(OBJRELDIR)/, $(REG_OBJS))

.PHONY: prep all clean debug release remake testdebug testrelease regdebug regrelease coverage profile

# Default build
all: debug release testdebug testrelease regdebug regrelease

#
# Debug rules
#
debug: $(LIBDBGEXE)

$(LIBDBGEXE): $(LIBDBGOBJS)
	$(CC) $(CFLAGS) $(SOFLAGS) $(DBGCFLAGS) -o $(LIBDBGEXE) $^ $(LDFLAGS)


#
# Release rules
#
release: $(LIBRELEXE)

$(LIBRELEXE): $(LIBRELOBJS)
	$(CC) $(CFLAGS) $(SOFLAGS) $(RELCFLAGS) -o $(LIBRELEXE) $^ $(LDFLAGS)

#
# Debug test rules
#
testdebug: debug $(TESTDBGEXE)

$(TESTDBGEXE): $(TESTDBGOBJS)
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(TESTDBGEXE) $^ $(DBGLDFLAGS)

#
# Release test rules
#
testrelease: release $(TESTRELEXE)

$(TESTRELEXE): $(TESTRELOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) -o $(TESTRELEXE) $^ $(RELLDFLAGS)

#
# Regression rules
#
regdebug: debug $(REGDBGEXE)

$(REGDBGEXE): $(REGDBGOBJS)
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(REGDBGEXE) $^ $(DBGLDFLAGS) $(LDFLAGS)
    
regrelease: release $(REGRELEXE)

$(REGRELEXE): $(REGRELOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) -o $(REGRELEXE) $^ $(RELLDFLAGS) $(LDFLAGS)

#
# Other rules
#
prep:
	@mkdir -p $(BINDIR) $(BINDBGDIR) $(BINRELDIR)
	@mkdir -p $(OBJDIR) $(OBJDBGDIR) $(OBJRELDIR)
	@mkdir -p $(LIBDIR) $(LIBDBGDIR) $(LIBRELDIR)
	@mkdir -p $(OBJLIBDIR) $(OBJLIBDBGDIR) $(OBJLIBRELDIR)

.ONESHELL:

coverage: regdebug
	@rm -f $(LIBDBGDIR)/*.gcda $(LIBDBGDIR)/*.gcda
	@LD_LIBRARY_PATH="$(LIBDBGDIR)" $(REGDBGEXE)
	@rm -f *.gcda *.gcno
	@cp $(OBJLIBDBGDIR)/*.gcno .
	@cp $(OBJLIBDBGDIR)/*.gcda .
	@gcov -n matrix.c
	@lcov -q -d . -c -o coverage.info > /dev/null
	@rm -rf $(COV_DIR)
	@genhtml -o $(COV_DIR) -t "couverture des tests" coverage.info > /dev/null
	@rm -f *.gc* coverage.info
	firefox $(COV_DIR)/index.html &

profile: regrelease
	mkdir -p $(PROF_DIR)
	rm -f $(PROF_DIR)/callgrind.out
	LD_LIBRARY_PATH="$(LIBRELDIR)" valgrind  --tool=callgrind --simulate-cache=yes --callgrind-out-file=$(PROF_DIR)/callgrind.out $(TESTRELEXE) 500
	kcachegrind $(PROF_DIR)/callgrind.out

clean:
	@rm -f  $(LIBRELEXE)  $(LIBRELOBJS)             \
	        $(LIBDBGEXE)  $(LIBDBGOBJS)             \
	        $(TESTDBGEXE) $(TESTDBGOBJS)            \
	        $(TESTRELEXE) $(TESTRELOBJS)            \
	        $(REGDBGEXE)  $(REGDBGOBJS)             \
	        $(REGRELEXE)  $(REGRELOBJS)             \
	        $(OBJLIBDBGDIR)/*.gc* $(OBJDBGDIR)/*.gc*
	@rm -fd $(OBJLIBRELDIR)  $(OBJLIBDBGDIR) $(OBJLIBDIR)   \
	        $(LIBRELDIR) $(LIBDBGDIR) $(LIBDIR)             \
	        $(OBJRELDIR) $(OBJDBGDIR) $(OBJDIR)             \
	        $(BINRELDIR) $(BINDBGDIR) $(BINDIR)                    
	@rm -rf $(COV_DIR)
	@rm -rf $(PROF_DIR)

