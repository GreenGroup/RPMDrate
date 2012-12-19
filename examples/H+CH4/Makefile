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
	$(F2PY) -c ch52008.f utility.f ch52008.pyf

pyf:
	$(F2PY) -m PES -h ch52008.pyf ch52008.f --overwrite-signature

clean:
	rm -f *.o PES.so
