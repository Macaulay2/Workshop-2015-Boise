bin_PROGRAMS = lrcalc schubmult

lrcalc_LDADD = liblrcalc.la
schubmult_LDADD = liblrcalc.la

lib_LTLIBRARIES = liblrcalc.la
liblrcalc_la_SOURCES = alloc.c hashtab.c list.c set.c vectarg.c vector.c \
                       lincomb.c maple.c symfcn.c schublib.c
liblrcalc_la_LDFLAGS = $(LIBLRCALC_LDFLAGS) -version-info 1:1:0

lrcalcincludedir = $(includedir)/lrcalc
lrcalcinclude_HEADERS = *.h


