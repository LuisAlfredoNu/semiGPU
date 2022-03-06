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
  CXX          := nvc++
endif

CXXFLAGS       += -std=c++11 -fPIC
OPTIMFLAGS      = -O2 -fast
OPENACCFLAGS := -acc -ta=$(DEVICE) -Minfo=accel

ARCHIVE      := ar
ARCHFLAG     := -rc

# OpenACC build flags
ifeq ($(OPENACC),1)
  CXXFLAGS     += $(OPENACCFLAGS)
  lLIBS        += -acc
  CXXFLAGS     += -DOPENACC_AVIL
endif

# Debug build flags
ifeq ($(DEBUGVERSION),1)
  CXXFLAGS     += -DDEBUG=1 -g -W -Wall -pedantic -Wno-long-long
  ifeq ($(OPENACC),1)
    CXXFLAGS   += -G -ta=$(DEVICE):lineinfo
    CXXFLAGS   += -DOPENACC_AVIL
  endif
else
  CXXFLAGS     += -DDEBUG=0 -Wall $(OPTIMFLAGS)
endif



INCDEFS        :=  #-include localdefs.h

# FILES
CPPOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
CCOBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
OBJS= $(CPPOBJS) $(CCOBJS)
HEADERS=$(wildcard *.h)
LIBOBJS=$(shell echo $(OBJS))
COMMONSTATICLIBNAME=$(shell cd ../..; pwd | sed -e 's;\(.*\)/\(.*\);\2;')
COMMONHEADERS=$(shell ls ../common/*.h)
COMMONSOURCES=$(shell ls ../common/*.cpp)
STATICLIB=$(TOP)/../common/lib$(COMMONSTATICLIBNAME).a
TESTEXECS=$(patsubst %.cc,%.x,$(wildcard *.cc))

# -L: FOLDER LIBRARY
	LFLAGS+=-L$(TOP)/../common/

# -l: LIBRARY
	lLIBS+=-l$(COMMONSTATICLIBNAME)
	lLIBS+=-pgf90libs
	lLIBS+=-llapack -lblas

# -I: INCLUDES
	IFLAGS+=-I. $(INCDEFS) -I$(TOP)/../common/

TARGET=target

$(TARGET): $(OBJS) $(STATICLIB) $(TESTEXECS)
	$(info LINKING: $@)

$(CPPOBJS): %.o: %.cpp
	$(info compiling: $@)
	@$(CXX) $(CXXFLAGS) -c $< -o $@  $(IFLAGS) 

$(CCOBJS): %.o: %.cc
	$(info compiling: $@)
	@$(CXX) $(CXXFLAGS) -c $<  -o $@  $(IFLAGS) 

$(TESTEXECS): %.x: %.o $(STATICLIB)
	$(info compiling: $@)
	$(CXX) $(LFLAGS) $< $(lLIBS) -o $@

$(STATICLIB): $(COMMONHEADERS) $(COMMONSOURCES)
	@cd $(TOP)/../common/;$(MAKE) DEBUGVERSION=$(DEBUGVERSION) OPENACC=$(OPENACC)

.PHONY: clean
clean:
	$(info CLEANING ALL)
	@$(RM) -f $(TARGET) $(OBJS) $(STATICLIB) $(TESTEXECS) 2>/dev/null || true

fullclean: clean
	$(info TOTAL CLEANING)
	@cd ../common/;make clean

libclean:
	 $(info CLEANING STATIC LIB)
	 @$(RM) -f  $(STATICLIB)

list:
	@echo $(OBJS)
	@echo $(HEADERS)

listobjs:
	@echo $(LIBOBJS)


