TOP      = $(shell pwd)

# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

# OS-specific build flags
ifneq ($(DARWIN),)
    CXXFLAGS   := -arch $(OS_ARCH)
else
  ifeq ($(OS_SIZE),32)
    CXXFLAGS   := -m32
  else
    CXXFLAGS   := -m64
  endif
endif

# COMPILERS OPTIONS
ifneq ($(DARWIN),)
  CXX          := g++-mp-4.9
else
  CXX          := g++
endif

CXXFLAGS       += -std=c++11 -fPIC
OPTIMFLAGS      = -O3 -funroll-loops -falign-loops=8 -ffast-math

ARCHIVE      := ar
ARCHFLAG     := -rc

# Debug build flags
ifeq ($(DEBUGVERSION),1)
  CXXFLAGS     += -DDEBUG=1 -g -W -Wall -pedantic -Wno-long-long
else
  CXXFLAGS     += -DDEBUG=0 -Wall $(OPTIMFLAGS)
endif

ifeq ($(USE_FLOATS_FOR_REALS),)
  CXXFLAGS     += -DUSE_FLOATS_FOR_REALS=0 
else
  CXXFLAGS     += -DUSE_FLOATS_FOR_REALS=$(USE_FLOATS_FOR_REALS)
endif

INCDEFS        := -include globaldefs.h

# FILES
CPPOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
OBJS= $(CPPOBJS)
HEADERS=$(wildcard *.h)
SOURCES=$(wildcard *.cpp)
LIBOBJS=$(shell echo $(OBJS))
COMMONSTATICLIBNAME=$(shell cd ../..; pwd | sed -e 's;\(.*\)/\(.*\);\2;')
STATICLIB="lib$(COMMONSTATICLIBNAME).a"

# -L: FOLDER LIBRARY
	LFLAGS+=#-L/usr/lib
	LFLAGS+=#-L.

# -l: LIBRARY
	lLIBS+=#-lgfnl_util
	lLIBS+= #-lapr-1 -laprutil-1 -llog4cxx
	lLIBS+=-llapack -lblas

# -I: INCLUDES
	IFLAGS+=-I. $(INCDEFS)

TARGET=target

$(TARGET): $(OBJS) $(HEADERS) $(SOURCES)
	$(info LINKING: $@)
	$(info ARCHIVING: $(LIBOBJS))
	@$(ARCHIVE) $(ARFLAGS) $(STATICLIB) $(LIBOBJS)

$(CPPOBJS): %.o: %.cpp
	$(info compiling: $@)
	@$(CXX) $(CXXFLAGS) -c $<  -o $@  $(IFLAGS) $(lLIBS)

.PHONY: clean
clean:
	$(info CLEANING ALL)
	@$(RM) -f $(TARGET) $(OBJS) $(STATICLIB) 2>/dev/null || true


list:
	@echo $(OBJS)
	@echo $(HEADERS)

listobjs:
	@echo $(LIBOBJS)

