ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src

lrcalcdatadir = $(datadir)/lrcalc
dist_lrcalcdata_DATA = lrcalc.maple README

TESTS = testsuite

lrcalc.maple: lrcalc.maple.src config.status
	cp lrcalc.maple.src lrcalc.maple
	echo "LRCALC_BIN_PATH := \"$(prefix)/bin/lrcalc\":" >> lrcalc.maple

