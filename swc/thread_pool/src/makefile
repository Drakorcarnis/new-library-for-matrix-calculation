#
# Coverage and profiling data
#
COV_DIR  = ../cov
PROF_DIR = ../prof
#
# Common compiler flags
#
CC      = gcc -pipe -fverbose-asm -g -ffast-math
AR      = ar rcs
CFLAGS  = -Wall -Werror -Wextra -std=c11
LDFLAGS = -lpthread
OPTFLAGS = -march=native -ffast-math

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
LIB_SRCS = thread_pool.c fifo.c
TEST_SRCS = test.c
REG_SRCS = regression.c
LIB_OBJS = $(LIB_SRCS:.c=.o)
TEST_OBJS = $(TEST_SRCS:.c=.o)
REG_OBJS = $(REG_SRCS:.c=.o)
TEST_EXE  = test
REGRESSION_EXE = regression
LIB  = thread_pool

#
# Lib Debug build settings
#
LIBDBGSHARED = $(LIBDBGDIR)/lib$(LIB).so
LIBDBGSTATIC = $(LIBDBGDIR)/lib$(LIB).a
LIBDBGOBJS = $(addprefix $(OBJLIBDBGDIR)/, $(LIB_OBJS))
DBGCFLAGS = -Og --coverage
DBGLDFLAGS = -L$(LIBDBGDIR) -l$(LIB) -lrt -lm
$(OBJLIBDBGDIR)/%.o: %.c .prep
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(DBGCFLAGS) -fPIC -o $@ $< -I$(INCLUDES)
$(OBJDBGDIR)/%.o: %.c .prep
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES)
#
# Lib Release build settings
#
LIBRELSHARED = $(LIBRELDIR)/lib$(LIB).so
LIBRELSTATIC = $(LIBRELDIR)/lib$(LIB).a
LIBRELOBJS = $(addprefix $(OBJLIBRELDIR)/, $(LIB_OBJS))
RELCFLAGS = -Ofast
RELLDFLAGS = -L$(LIBRELDIR) -l$(LIB) -lrt -lm -lgomp
$(OBJLIBRELDIR)/%.o: %.c .prep
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(RELCFLAGS) -fPIC -o $@ $< -I$(INCLUDES) 
$(OBJRELDIR)/%.o: %.c .prep
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(RELCFLAGS) -o $@ $< -I$(INCLUDES)
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
	@touch .prep

debug: $(LIBDBGSHARED) $(LIBDBGSTATIC)

$(LIBDBGSHARED): $(LIBDBGOBJS)
	$(CC) $(CFLAGS) $(OPTFLAGS) $(DBGCFLAGS) -shared -o $(LIBDBGSHARED) $^ $(LDFLAGS)

$(LIBDBGSTATIC): $(LIBDBGOBJS)
	$(AR) -o $(LIBDBGSTATIC) $(LIBDBGOBJS)

release: $(LIBRELSHARED) $(LIBRELSTATIC)

$(LIBRELSHARED): $(LIBRELOBJS)
	$(CC) $(CFLAGS) $(OPTFLAGS) $(RELCFLAGS) -shared -o $(LIBRELSHARED) $^ $(LDFLAGS)

$(LIBRELSTATIC): $(LIBRELOBJS)
	$(AR) -o $(LIBRELSTATIC) $(LIBRELOBJS)

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

#
# Other rules
#

.ONESHELL:

coverage: regdebug
	@rm -f $(LIBDBGDIR)/*.gcda $(LIBDBGDIR)/*.gcda
	@LD_LIBRARY_PATH="$(LIBDBGDIR)" $(REGDBGEXE)
	@rm -f *.gcda *.gcno
	@cp $(OBJLIBDBGDIR)/*.gcno .
	@cp $(OBJLIBDBGDIR)/*.gcda .
	@gcov -n thread_pool.c
	@lcov -q -d . -c -o coverage.info > /dev/null
	@rm -rf $(COV_DIR)
	@genhtml -o $(COV_DIR) -t "couverture des tests" coverage.info > /dev/null
	@rm -f *.gc* coverage.info
	firefox $(COV_DIR)/index.html &

profile: regrelease
	mkdir -p $(PROF_DIR)
	rm -f $(PROF_DIR)/callgrind.out
	LD_LIBRARY_PATH="$(LIBRELDIR)" valgrind  --tool=callgrind --dump-instr=yes --simulate-cache=yes --callgrind-out-file=$(PROF_DIR)/callgrind.out $(TESTRELEXE) 20000 20
	kcachegrind $(PROF_DIR)/callgrind.out

clean:
	@rm -f  .prep
	@rm -f  $(LIBRELSHARED)  $(LIBDBGSHARED)        \
	        $(LIBRELSTATIC)  $(LIBDBGSTATIC)        \
	        $(LIBRELOBJS)  $(LIBDBGOBJS)            \
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

