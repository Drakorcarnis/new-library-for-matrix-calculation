#
# Compiler flags
#
CC     = gcc
CFLAGS = -Wall -Werror -Wextra
LDFLAGS = -lm
TESTDBGLDFLAGS = -L$(DBGDIR) -lmatrix
TESTRELLDFLAGS = -L$(RELDIR) -lmatrix

#
# Project files
#
SRCS = matrix.c
TEST_SRCS = test.c
OBJS = $(SRCS:.c=.o)
TEST_OBJS = $(TEST_SRCS:.c=.o)
EXE  = test
SO  = libmatrix.so

#
# Debug build settings
#
DBGDIR = ../bin/debug
DBGEXE = $(DBGDIR)/$(SO)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGCFLAGS = -shared -fPIC -g -O0 -DDEBUG

#
# Release build settings
#
RELDIR = ../bin/release
RELEXE = $(RELDIR)/$(SO)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELCFLAGS = -shared -fPIC -O3 -DNDEBUG

#
# Test Debug build settings
#
TESTDBGDIR = ../test/debug
TESTDBGEXE = $(TESTDBGDIR)/$(EXE)
TESTDBGOBJS = $(addprefix $(TESTDBGDIR)/, $(TEST_OBJS))
TESTDBGCFLAGS = -g -O0 -DDEBUG

#
# Test Release build settings
#
TESTRELDIR = ../test/release
TESTRELEXE = $(TESTRELDIR)/$(EXE)
TESTRELOBJS = $(addprefix $(TESTRELDIR)/, $(TEST_OBJS))
TESTRELCFLAGS = -O3 -DNDEBUG

.PHONY: all clean debug prep release remake testdebug testrelease

# Default build
all: prep release debug testdebug testrelease

#
# Debug rules
#
debug: $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(DBGEXE) $^ $(LDFLAGS)

$(DBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $<

#
# Release rules
#
release: $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) -o $(RELEXE) $^ $(LDFLAGS)

$(RELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(RELCFLAGS) -o $@ $<
#
# Debug test rules
#
testdebug: $(TESTDBGEXE)

$(TESTDBGEXE): $(TESTDBGOBJS)
	$(CC) $(CFLAGS) $(TESTDBGCFLAGS) -o $(TESTDBGEXE) $^ $(TESTDBGLDFLAGS)
$(TESTDBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(TESTDBGCFLAGS) -o $@ $<

#
# Release test rules
#
testrelease: $(TESTRELEXE)

$(TESTRELEXE): $(TESTRELOBJS)
	$(CC) $(CFLAGS) $(TESTRELCFLAGS) -o $(TESTRELEXE) $^ $(TESTRELLDFLAGS)

$(TESTRELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(TESTRELCFLAGS) -o $@ $<

#
# Other rules
#
prep:
	@mkdir -p $(DBGDIR) $(RELDIR) $(TESTDBGDIR) $(TESTRELDIR)

remake: clean all

clean:
	rm -f $(RELEXE) $(RELOBJS) $(DBGEXE) $(DBGOBJS) $(TESTDBGEXE) $(TESTDBGOBJS) $(TESTRELEXE) $(TESTRELOBJS)