#
# Common compiler flags
#
CC     = gcc
CFLAGS = -Wall -Werror -Wextra
LDFLAGS = -lm
SOFLAGS = -shared -fPIC
BINDIR = ../bin

#
# Project files
#
INCLUDES = includes
SRCS = matrix.c matrix_tools.c
TEST_SRCS = test.c
REG_SRCS = regression.c
OBJS = $(SRCS:.c=.o)
TEST_OBJS = $(TEST_SRCS:.c=.o)
REG_OBJS = $(REG_SRCS:.c=.o)
TEST_EXE  = test
REGRESSION_EXE = regression
SO  = libmatrix.so

#
# Debug build settings
#
DBGDIR = $(BINDIR)/debug
LIBDBGDIR = $(DBGDIR)/lib
LIBDBGEXE = $(LIBDBGDIR)/$(SO)
LIBDBGOBJS = $(addprefix $(LIBDBGDIR)/, $(OBJS))
DBGCFLAGS = -g -O0 -DDEBUG --coverage
DBGLDFLAGS = -L$(DBGDIR)/lib -lmatrix
$(LIBDBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(SOFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES)
$(DBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $< -I$(INCLUDES)
#
# Release build settings
#
RELDIR = $(BINDIR)/release
LIBRELDIR = $(RELDIR)/lib
LIBRELEXE = $(LIBRELDIR)/$(SO)
LIBRELOBJS = $(addprefix $(LIBRELDIR)/, $(OBJS))
RELCFLAGS = -O3 -DNDEBUG
RELLDFLAGS = -L$(RELDIR)/lib -lmatrix
$(LIBRELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(SOFLAGS) $(RELCFLAGS) -o $@ $< -I$(INCLUDES)
$(RELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(RELCFLAGS) -o $@ $< -I$(INCLUDES)

#
# Test Debug build settings
#
TESTDBGEXE = $(DBGDIR)/$(TEST_EXE)
TESTDBGOBJS = $(addprefix $(DBGDIR)/, $(TEST_OBJS))

#
# Test Release build settings
#
TESTRELEXE = $(RELDIR)/$(TEST_EXE)
TESTRELOBJS = $(addprefix $(RELDIR)/, $(TEST_OBJS))

#
# Regression build settings
#
REGEXE = $(DBGDIR)/$(REGRESSION_EXE)
REGOBJS = $(addprefix $(DBGDIR)/, $(REG_OBJS))

.PHONY: all clean debug prep release remake testdebug testrelease reg

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
testdebug: $(TESTDBGEXE)

$(TESTDBGEXE): $(TESTDBGOBJS) debug
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(TESTDBGEXE) $< $(DBGLDFLAGS)

#
# Release test rules
#
testrelease: $(TESTRELEXE)

$(TESTRELEXE): $(TESTRELOBJS) release
	$(CC) $(CFLAGS) $(RELCFLAGS) -o $(TESTRELEXE) $< $(RELLDFLAGS)

#
# Regression rules
#
reg: $(REGEXE)

$(REGEXE): $(REGOBJS) debug
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(REGEXE) $< $(DBGLDFLAGS) $(LDFLAGS)

#
# Other rules
#
prep:
	@mkdir -p $(DBGDIR) $(RELDIR) $(LIBDBGDIR) $(LIBRELDIR)

remake: clean all

lcov: reg
	@LD_LIBRARY_PATH="$(LIBDBGDIR)" $(REGEXE)
	@rm -f *.gcda
	@rm -f *.gcno
	@cp $(LIBDBGDIR)/*.gc* .
	@cp $(DBGDIR)/*.gc* .
	@gcov -n matrix.c
	@lcov -q -d . -c -o coverage.info > /dev/null
	@genhtml -o ../cov -t "couverture des tests" coverage.info > /dev/null
	@rm -f *.gc* coverage.info
	firefox ../cov/index.html

clean:
	rm -f $(LIBRELEXE)  $(LIBRELOBJS)                       \
	      $(LIBDBGEXE)  $(LIBDBGOBJS)  $(LIBDBGDIR)/*.gc*   \
	      $(TESTDBGEXE) $(TESTDBGOBJS) $(TESTDBGEXE).gc*    \
	      $(TESTRELEXE) $(TESTRELOBJS)                      \
	      $(REGEXE)     $(REGOBJS)     $(REGEXE).gc*  
	rm -fd $(LIBDBGDIR)  $(LIBRELDIR)                        \
	       $(DBGDIR) $(RELDIR)                               \
	       $(BINDIR)
	rm -rf ../cov

