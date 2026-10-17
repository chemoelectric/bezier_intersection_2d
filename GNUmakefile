.DELETE_ON_ERROR:
.DEFAULT_TARGET := default

builddir = ___build___

# GNU m4 or OpenBSD m4. Those are the only two I know how to or wish
# to support. OpenBSD is surely more secure, thanks to the magic of
# simplicity, but GNU is far, far, far, far, far more capable, thanks
# to recursion optimization. (That is, if you want to prove you can
# write <fill in the blank> in m4, use GNU m4. I was able to do the
# Rosetta Code “Evolutionary Algorithm” task in GNU m4, though it
# required a long run time on the AMD Opterons I had then. I would
# probably prefer OpenBSD stuff, were I a network admin, but I am
# instead a microcomputer hobbyist who started on a TRS-89 in 1979. I
# like GNU for my own stuff. One never fills an arbitrarily sized
# array. There needs to be a more basic way of doing this. Everything
# is wrong in the world of computer programming.)
GNU_M4 = m4

M4 = $(GNU_M4)
M4FLAGS = -P
TESTM4FLAGS = -s -P

CC = gcc
CFLAGS = -O3 -g -std=gnu99 -I$(builddir)
LIBS = -lm

UNCRUSTIFY = uncrustify
UNCRUSTIFYFLAGS = -q -c uncrustify.cfg --replace

.PHONY: default all clean veryclean
default:: all

all:: bezier_intersection_2d.h
veryclean:: ; -rm -f bezier_intersection_2d.h

%.h: %.h.m4
	$(M4) $(M4FLAGS) $(<) > $(@)
	$(UNCRUSTIFY) $(UNCRUSTIFYFLAGS) $(@)

$(builddir)/test_%.h: %.h.m4 | $(builddir)
	$(M4) $(TESTM4FLAGS) $(<) > $(@)

$(builddir): GNUmakefile
	mkdir -p $(@)
clean:: ; -rm -R -f $(builddir)/*.{c,h,o}

.PHONY: indent_test
all:: indent_test
indent_test:
	$(UNCRUSTIFY) $(UNCRUSTIFYFLAGS) test_bezier_intersection_2d.c

$(builddir)/test_bezier_intersection_2d: \
		test_bezier_intersection_2d.c \
		$(builddir)/test_bezier_intersection_2d.h
	$(CC) $(CFLAGS) $(<) $(LIBS) -o $(@)
clean:: ; -rm -f $(builddir)/test_bezier_intersection_2d

.PHONY: test check
test check:: $(builddir)/test_bezier_intersection_2d
	$(builddir)/test_bezier_intersection_2d
