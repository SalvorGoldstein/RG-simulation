CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wno-unused-parameter -O2 \
           -MMD -MP \
           -I . \
           -I rgsim/third_party/autodiff \
           -I rgsim/third_party/eigen-3.4.0

SOURCES = main.cc example_metric.cc
OBJECTS = $(SOURCES:.cc=.o)
DEPS    = $(OBJECTS:.o=.d)

TARGET = main

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

execute: $(TARGET)
	mkdir -p data
	./$(TARGET)

aff:
	python python/Geodesics_topview.py

clean:
	rm -f *.o *.d $(TARGET)

-include $(DEPS)

.PHONY: all execute aff clean