.ONESHELL:
#
# Coverage and profiling data
#

COV_DIR = ../cov
PROF_DIR = ../prof
#
# Common compiler flags
#
CC      = gcc -pipe -fverbose-asm -g -ffast-math
CFLAGS  = -Wall -Werror -Wextra -std=c11
LDFLAGS = -lm -lgomp
SOFLAGS = -shared -fPIC -fopenmp -march=native -ffast-math

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
LIB_SRCS = matrix.c tools.c plu.c cholesky.c check.c raw.c
TEST_SRCS = test.c
REG_SRCS = regression.c
LIB_OBJS = $(LIB_SRCS:.c=.o)
TEST_OBJS = $(TEST_SRCS:.c=.o)
REG_OBJS = $(REG_SRCS:.c=.o)
TEST_EXE  = test
REGRESSION_EXE = regression
LIB  = matrix

#
# SWC thread_pool lib build settings
#
THPOOL_DIR = ../swc/thread_pool
THPOOL_INCLUDES = $(THPOOL_DIR)/src/includes
THPOOL_LIBDIR = $(THPOOL_DIR)/lib
THPOOL_LIBDBGDIR = $(THPOOL_LIBDIR)/debug
THPOOL_LIBRELDIR = $(THPOOL_LIBDIR)/release
THPOOL_LIB = thread_pool
THPOOL_LIB_DBG_STATIC = $(THPOOL_LIBDBGDIR)/lib$(THPOOL_LIB).a
THPOOL_LIB_REL_STATIC = $(THPOOL_LIBRELDIR)/lib$(THPOOL_LIB).a
#
# Lib Debug build settings
#
LIBDBGSHARED = $(LIBDBGDIR)/lib$(LIB).so
LIBDBGOBJS = $(addprefix $(OBJLIBDBGDIR)/, $(LIB_OBJS))
DBGCFLAGS = -Og --coverage 
DBGLDFLAGS = -L$(LIBDBGDIR) -l$(LIB)
$(OBJLIBDBGDIR)/%.o: %.c .prep
	$(CC) -c $(CFLAGS) $(SOFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES) -I$(THPOOL_INCLUDES)
$(OBJDBGDIR)/%.o: %.c .prep
	$(CC) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES)
#
# Lib Release build settings
#
LIBRELSHARED = $(LIBRELDIR)/lib$(LIB).so
LIBRELOBJS = $(addprefix $(OBJLIBRELDIR)/, $(LIB_OBJS))
RELCFLAGS = -Ofast -fopt-info-vec-optimized
RELLDFLAGS = -L$(LIBRELDIR) -l$(LIB)
$(OBJLIBRELDIR)/%.o: %.c .prep
	$(CC) -c $(CFLAGS) $(SOFLAGS) $(RELCFLAGS) -o $@ $< -I$(INCLUDES) -I$(THPOOL_INCLUDES)
$(OBJRELDIR)/%.o: %.c .prep
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

.PHONY: all clean debug release remake testdebug testrelease regdebug regrelease coverage profile

# Default build
all: debug release testdebug testrelease regdebug regrelease

.prep:
	@mkdir -p $(BINDIR) $(BINDBGDIR) $(BINRELDIR)
	@mkdir -p $(OBJDIR) $(OBJDBGDIR) $(OBJRELDIR)
	@mkdir -p $(LIBDIR) $(LIBDBGDIR) $(LIBRELDIR)
	@mkdir -p $(OBJLIBDIR) $(OBJLIBDBGDIR) $(OBJLIBRELDIR)
	touch .prep

$(THPOOL_LIB_REL_STATIC):
	@(cd $(THPOOL_DIR)/src && $(MAKE) release)

$(THPOOL_LIB_DBG_STATIC):
	@(cd $(THPOOL_DIR)/src && $(MAKE) debug)
    
debug: $(LIBDBGSHARED)

$(LIBDBGSHARED): $(THPOOL_LIB_DBG_STATIC) $(LIBDBGOBJS)
	$(CC) $(CFLAGS) $(SOFLAGS) $(DBGCFLAGS) -o $(LIBDBGSHARED) $(LIBDBGOBJS) $(THPOOL_LIB_DBG_STATIC) $(LDFLAGS)

release: $(LIBRELSHARED)

$(LIBRELSHARED): $(THPOOL_LIB_REL_STATIC) $(LIBRELOBJS)
	$(CC) $(CFLAGS) $(SOFLAGS) $(RELCFLAGS) -o $(LIBRELSHARED) $(LIBRELOBJS)  $(THPOOL_LIB_REL_STATIC) $(LDFLAGS)

testdebug: $(TESTDBGEXE)

$(TESTDBGEXE): $(LIBDBGSHARED) $(TESTDBGOBJS)
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(TESTDBGEXE) $(TESTDBGOBJS) $(DBGLDFLAGS)

testrelease: $(TESTRELEXE)

$(TESTRELEXE): $(LIBRELSHARED) $(TESTRELOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) -o $(TESTRELEXE) $(TESTRELOBJS) $(RELLDFLAGS)

regdebug: $(REGDBGEXE)

$(REGDBGEXE): $(LIBDBGSHARED) $(REGDBGOBJS)
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(REGDBGEXE) $(REGDBGOBJS) $(DBGLDFLAGS) $(LDFLAGS)
  
regrelease: $(REGRELEXE)

$(REGRELEXE): $(LIBRELSHARED) $(REGRELOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) -o $(REGRELEXE) $(REGRELOBJS) $(RELLDFLAGS) $(LDFLAGS)


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
	LD_LIBRARY_PATH="$(LIBRELDIR)" valgrind  --tool=callgrind --simulate-cache=yes --callgrind-out-file=$(PROF_DIR)/callgrind.out $(TESTRELEXE) 1000
	kcachegrind $(PROF_DIR)/callgrind.out

clean:
	@rm -f  .prep
	@(cd $(THPOOL_DIR)/src && $(MAKE) $@)
	@rm -f  $(LIBRELSHARED)  $(LIBRELOBJS)             \
	        $(LIBDBGSHARED)  $(LIBDBGOBJS)             \
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

