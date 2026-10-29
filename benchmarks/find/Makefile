# homer2 only, with the flags and link line of HOMER v5.1's cpp/Makefile.
#	make -f benchmarks/find/Makefile -C <HOMER tree> OUT=<dir>    builds <dir>/homer2

CXX = g++
DEBUG = -O3 -Wall -g
LIBRARIES = -lpthread -lstdc++ -lm
OUT ?= build
SRC = cpp

OBJS = $(OUT)/homer2.o $(OUT)/Motif2.o $(OUT)/Hashtable.o $(OUT)/statistics.o
HEADERS = $(wildcard $(SRC)/*.h)

$(OUT)/homer2 : $(OBJS)
	$(CXX) $(OBJS) $(LIBRARIES) $(DEBUG) -o $@

$(OUT)/%.o : $(SRC)/%.cpp $(HEADERS) | $(OUT)
	$(CXX) $< -c $(DEBUG) -o $@

$(OUT) :
	mkdir -p $(OUT)
