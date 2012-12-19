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
	$(F2PY) -c c2h7.f c2h7.pyf

pyf:
	$(F2PY) -m PES -h c2h7.pyf c2h7.f --overwrite-signature

clean:
	rm -f *.o *.so
