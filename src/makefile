#
# Coverage and profiling data
#
COV_DIR = ../cov
PROF_DIR = ../prof
#
# Common compiler flags
#
CC     = gcc
CFLAGS = -Wall -Werror -Wextra
LDFLAGS = -lm -lgomp
SOFLAGS = -shared -fPIC -fopenmp -mavx

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
LIB_SRCS = matrix.c matrix_tools.c
TEST_SRCS = test.c
REG_SRCS = regression.c
LIB_OBJS = $(LIB_SRCS:.c=.o)
TEST_OBJS = $(TEST_SRCS:.c=.o)
REG_OBJS = $(REG_SRCS:.c=.o)
TEST_EXE  = test
REGRESSION_EXE = regression
LIB_SO  = libmatrix.so

#
# Lib Debug build settings
#
LIBDBGEXE = $(LIBDBGDIR)/$(LIB_SO)
LIBDBGOBJS = $(addprefix $(OBJLIBDBGDIR)/, $(LIB_OBJS))
DBGCFLAGS = -g -O0 -DDEBUG --coverage
DBGLDFLAGS = -L$(LIBDBGDIR) -lmatrix
$(OBJLIBDBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(SOFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES)
$(OBJDBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES)
#
# Lib Release build settings
#
LIBRELEXE = $(LIBRELDIR)/$(LIB_SO)
LIBRELOBJS = $(addprefix $(OBJLIBRELDIR)/, $(LIB_OBJS))
RELCFLAGS = -O3 -DNDEBUG
RELLDFLAGS = -L$(LIBRELDIR) -lmatrix
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
REGEXE = $(BINDBGDIR)/$(REGRESSION_EXE)
REGOBJS = $(addprefix $(OBJDBGDIR)/, $(REG_OBJS))

.PHONY: all clean debug prep release remake testdebug testrelease reg coverage profile

# Default build
all: prep release debug testdebug testrelease reg

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
reg: debug $(REGEXE)

$(REGEXE): $(REGOBJS)
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(REGEXE) $^ $(DBGLDFLAGS) $(LDFLAGS)

#
# Other rules
#
prep:
	@mkdir -p $(BINDIR) $(BINDBGDIR) $(BINRELDIR)
	@mkdir -p $(OBJDIR) $(OBJDBGDIR) $(OBJRELDIR)
	@mkdir -p $(LIBDIR) $(LIBDBGDIR) $(LIBRELDIR)
	@mkdir -p $(OBJLIBDIR) $(OBJLIBDBGDIR) $(OBJLIBRELDIR)

remake: clean all

.ONESHELL:

coverage: reg
	@rm -f $(LIBDBGDIR)/*.gcda $(LIBDBGDIR)/*.gcda
	@LD_LIBRARY_PATH="$(LIBDBGDIR)" $(REGEXE)
	@rm -f *.gcda *.gcno
	@cp $(OBJLIBDBGDIR)/*.gcno .
	@cp $(OBJLIBDBGDIR)/*.gcda .
	@gcov -n matrix.c
	@lcov -q -d . -c -o coverage.info > /dev/null
	@rm -rf $(COV_DIR)
	@genhtml -o $(COV_DIR) -t "couverture des tests" coverage.info > /dev/null
	@rm -f *.gc* coverage.info
	firefox $(COV_DIR)/index.html &

profile: reg
	mkdir -p $(PROF_DIR)
	rm -f $(PROF_DIR)/callgrind.out
	LD_LIBRARY_PATH="$(LIBDBGDIR)" valgrind  --tool=callgrind --simulate-cache=yes --callgrind-out-file=$(PROF_DIR)/callgrind.out $(REGEXE)
	kcachegrind $(PROF_DIR)/callgrind.out

clean:
	@rm -f  $(LIBRELEXE)  $(LIBRELOBJS)             \
	        $(LIBDBGEXE)  $(LIBDBGOBJS)             \
	        $(TESTDBGEXE) $(TESTDBGOBJS)            \
	        $(TESTRELEXE) $(TESTRELOBJS)            \
	        $(REGEXE)     $(REGOBJS)                \
	        $(OBJLIBDBGDIR)/*.gc* $(OBJDBGDIR)/*.gc*
	@rm -fd $(OBJLIBRELDIR)  $(OBJLIBDBGDIR) $(OBJLIBDIR)   \
	        $(LIBRELDIR) $(LIBDBGDIR) $(LIBDIR)             \
	        $(OBJRELDIR) $(OBJDBGDIR) $(OBJDIR)             \
	        $(BINRELDIR) $(BINDBGDIR) $(BINDIR)                    
	@rm -rf $(COV_DIR)
	@rm -rf $(PROF_DIR)

