CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wno-unused-parameter -O2 -I cc -I cc/third_party/autodiff -I cc/third_party/eigen-3.4.0

SOURCES1 = cc/multiplegeo.cc cc/geodesic.cc cc/funcoperator.cc cc/metricoperator.cc cc/metric.cc
SOURCES2 = cc/forvideo.cc cc/geodesic.cc cc/funcoperator.cc cc/metricoperator.cc cc/metric.cc
OBJETS1  = $(SOURCES1:.cc=.o)
OBJETS2  = $(SOURCES2:.cc=.o)
HEADERS = cc/types.h cc/geodesic.h cc/funcoperator.h cc/metricoperator.h cc/metric.h



multiplegeo: $(OBJETS1)
	$(CXX) $(CXXFLAGS) $(OBJETS1) -o multiplegeo

: $(OBJETS2)
	$(CXX) $(CXXFLAGS) $(OBJETS1) -o forvideo

%.o: %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

execute1: multiplegeo
	./multiplegeo

execute2: forvideo
	./forvideo

aff1: Geodesics.dat
	python python/Schwartzschield_aff.py

aff2: data
	python python/Geodesics_topview.py

Geodesics1.dat: multiplegeo
	./multiplegeo

clean:
	rm -f $(OBJETS1) multiplegeo

.PHONY: all execute aff clean