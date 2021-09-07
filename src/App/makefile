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
  CXXFLAGS     += -DDEBUG=0 -Wall -pedantic $(OPTIMFLAGS)
endif

INCDEFS        :=  -include ../common/globaldefs.h

# FILES
CPPOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
CCOBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
OBJS= $(CPPOBJS) $(CCOBJS)
HEADERS=$(wildcard *.h)
SOURCES=$(wildcard *.cpp *.cc)
LIBOBJS=$(shell echo $(OBJS))
COMMONSTATICLIBNAME=$(shell cd ../..; pwd | sed -e 's;\(.*\)/\(.*\);\2;')
COMMONHEADERS=$(shell ls ../common/*.h)
COMMONSOURCES=$(shell ls ../common/*.cpp)
STATICLIB=$(TOP)/../common/lib$(COMMONSTATICLIBNAME).a
TESTEXECS=$(patsubst %.cc,%.x,$(wildcard *.cc))

# -L: FOLDER LIBRARY
	LFLAGS+=-L$(TOP)/../common/
	LFLAGS+=#-L.

# -l: LIBRARY
	lLIBS+=-l$(COMMONSTATICLIBNAME)
	lLIBS+=-llapack -lblas

# -I: INCLUDES
	IFLAGS+=-I. $(INCDEFS) -I$(TOP)/../common/

TARGET=target

$(TARGET): $(OBJS) $(TESTEXECS) $(STATICLIB)
	$(info LINKING: $@)

$(CPPOBJS): %.o: %.cpp
	$(info compiling: $@)
	@$(CXX) $(CXXFLAGS) -c $< -o $@  $(IFLAGS) 

$(CCOBJS): %.o: %.cc
	$(info compiling: $@)
	@$(CXX) $(CXXFLAGS) -c $<  -o $@  $(IFLAGS) 

$(TESTEXECS): %.x: %.o $(STATICLIB) $(CPPOBJS) $(SOURCES) $(HEADERS) $(COMMONHEADERS) $(COMMONSOURCES)
	$(info compiling: $@)
	$(CXX) $(LFLAGS) $< $(CPPOBJS) $(lLIBS) -o $@

$(STATICLIB): $(COMMONHEADERS) $(COMMONSOURCES)
	@cd $(TOP)/../common/;make

.PHONY: clean
clean:
	$(info CLEANING ALL)
	@$(RM) -f $(TARGET) $(OBJS) $(STATICLIB) $(TESTEXECS) 2>/dev/null || true

fullclean: clean
	@cd ../common/;make clean


list:
	@echo $(OBJS)
	@echo $(HEADERS)

listobjs:
	@echo $(LIBOBJS)

