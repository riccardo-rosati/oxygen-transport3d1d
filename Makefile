# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2016-2017
#
#                    Copyright D. Brambilla 2016
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for library installation
#   AUTHOR      : Stefano Brambilla <s.brambilla93@gmail.com>
#   DATE        : September 2016
# ====================================================================
include ../config.mk

CPPFLAGS+=-I. -I$(GETFEM_PREFIX)/include -I$(mkBoostInc) -I$(PROBLEM_FLUID_PREFIX)/include
CXXFLAGS+=-std=c++11 -Wall -fPIC


ifeq ($(DEBUG),yes)
  OPTFLAGS=-g -Wall
else
  OPTFLAGS=-O3 -march=native
  CPPFLAGS+=-DNDEBUG
endif

ifeq ($(VERBOSE),yes)
  CPPFLAGS+=-DM3D1D_VERBOSE_
endif

LIBDIR=../lib
LIBNAME=transport3d1d
LIBFILE=lib$(LIBNAME).so
LIBSTATIC=lib$(LIBNAME).a

LIBSRC=$(wildcard *.cpp)
LIBOBJS=$(LIBSRC:.cpp=.o)
LIBHEADERS=$(wildcard *.hpp)

.PHONY: all clean distclean library

all: library
	@echo
	@echo Library installed!

library: $(LIBOBJS)
	$(CXX) -shared -o $(LIBDIR)/$(LIBFILE) $(LIBOBJS)

static: $(LIBOBJS)
	install -d $(LIBDIR)
	ar -r $(LIBDIR)/$(LIBSTATIC) $(LIBOBJS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -o $@ -c $<

clean:
	$(RM) $(LIBOBJS) *~

distclean: clean
	$(RM) $(LIBDIR)/$(LIBFILE) $(LIBDIR)/$(LIBSTATIC) 
