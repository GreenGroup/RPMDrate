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
	$(F2PY) -c ch4oh.f utility.f ch4oh.pyf

pyf:
	$(F2PY) -m PES -h ch4oh.pyf ch4oh.f --overwrite-signature

clean:
	rm -f *.o *.so
