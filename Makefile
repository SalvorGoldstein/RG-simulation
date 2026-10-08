CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wno-unused-parameter -O2 -I cc -I cc/third_party/autodiff -I cc/third_party/eigen-3.4.0

SOURCES = cc/main.cc cc/geodesic.cc cc/funcoperator.cc cc/metricoperator.cc cc/metric.cc
OBJETS  = $(SOURCES:.cc=.o)
HEADERS = cc/types.h cc/geodesic.h cc/funcoperator.h cc/metricoperator.h cc/metric.h

all: main

main: $(OBJETS)
	$(CXX) $(CXXFLAGS) $(OBJETS) -o main

%.o: %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

execute: main
	./main

aff: Geodesics.dat
	python python/Schwartzschield_aff.py

Geodesics.dat: main
	./main

clean:
	rm -f $(OBJETS) main

.PHONY: all execute aff clean