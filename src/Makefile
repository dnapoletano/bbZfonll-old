CXXFLAGS =  $(shell lhapdf-config --cflags)
CXXFLAGS += $(shell gsl-config --cflags)
LDFLAGS  = $(shell lhapdf-config --libs)
LDFLAGS  += $(shell gsl-config --libs)  

LDFLAGS += -L. -lbbh -L$(WORK)/lib -lgfortran  #-L/usr/lib64 -lgfortran

PYPREFIX= $(shell python3-config --prefix)
PYFLAGS = ${CXXFLAGS}
PYFLAGS += $(shell python3-config --cflags)
PYFLAGS += $(shell python3-config --ldflags)

PYFLAGS  += ${LDFLAGS}

INSTALLPATH=plots_HXSWG

CXX = g++

all: libbbh.a test

libbbh.a: bits.o # bbh.o bbhnnloexact.o
	ar r $@ $^

bits.o: bbznnlo/bits.f # bbh.o bbhnnloexact.o: bbhnnlo/bbh.f # bbhnnlo/bbhnnloexact.f
	gfortran -c -g -fPIC $^

%.o: %.C
	$(CXX) -fPIC -c -g $< ${CXXFLAGS} 

test:  4F.o Calculate.o 4F0_LO.o 4F0_NLO.o 5F.o Massless_lim.o PDF.o 
	$(CXX) -g -Wall -Werror -shared $^ ${LDFLAGS} -o libtest.so

bbHFONLL_wrap.cxx: bbHFONLL.i
	swig -c++ -python $<

python: 4F.C Calculate.C 4F0_LO.C 4F0_NLO.C 5F.C Massless_lim.C PDF.C bbHFONLL_wrap.cxx
	$(CXX) -g -fPIC -o _bbHFONLL.so $^ ${PYFLAGS} -Wl,-rpath,${PYPREFIX}/lib

install: _bbHFONLL.so bbHFONLL.py
	cp $^ ${INSTALLPATH}/.

clean:
	rm -rf *.o test *~ ../lib/libbbh.a *_wrap.cxx ../lib/_*.so ../lib/bbHFONLL.py *.dSYM libbbh.a _*.so bbHFONLL.py

