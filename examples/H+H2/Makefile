################################################################################
#
#   This Makefile will compile the potential energy surface module to expose it
#   to Python using f2py.
#
################################################################################

F2PY=f2py

.PHONY: all build pyf clean

all: build

build:
	$(F2PY) -c h3_bkmp2.f h3_bkmp2.pyf

pyf:
	$(F2PY) -m PES -h h3_bkmp2.pyf h3_bkmp2.f --overwrite-signature

clean:
	rm -f *.o *.so
