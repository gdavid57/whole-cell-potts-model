# Makefile optimized for maximum performance
CXX = g++

# Aggressive optimization flags
CXXFLAGS_PERFORMANCE = -m64 -O2 -DNDEBUG -march=native -mtune=native \
                      -ffast-math -funroll-loops -finline-functions \
                      -fomit-frame-pointer -Wall -Wextra

# Debug flags (original)
CXXFLAGS_DEBUG = -m64 -fdiagnostics-color=always -g -fsanitize=address -fsanitize=undefined -Wall

# Libraries
# LDLIBS = -pthread

# Source and target
SRC = main.cpp cell_state.cpp
TARGET = wcm

# Python binding configuration
PYTHON_BINDING_SRC = py_cell_state.cpp cell_state.cpp
PYTHON_MODULE = py_cell_state$(shell python3-config --extension-suffix)
PYTHON_FLAGS = $(shell python3 -m pybind11 --includes) $(shell python3-config --ldflags)

# Rules
.PHONY: all performance debug clean python

# Default build (performance an python)
all: performance python

# Maximum performance build
performance:
	$(CXX) $(CXXFLAGS_PERFORMANCE) $(SRC) -o $(TARGET)
	@echo "=== OPTIMIZED PERFORMANCE BUILD ==="
	@echo "Flags used: $(CXXFLAGS_PERFORMANCE)"

# Debug build (original)
debug:
	$(CXX) $(CXXFLAGS_DEBUG) $(SRC) -o $(TARGET)

# Python binding build
python:
	$(CXX) $(CXXFLAGS_PERFORMANCE) $(PYTHON_FLAGS) -shared -fPIC $(PYTHON_BINDING_SRC) -o $(PYTHON_MODULE)
	@echo "=== PYTHON BINDING BUILD ==="
	@echo "Module created: $(PYTHON_MODULE)"
	@echo "Flags used: $(CXXFLAGS_PERFORMANCE) $(PYTHON_FLAGS)"

# Cleanup
clean:
	rm -f $(TARGET) *.o $(PYTHON_MODULE)

# Quick performance test
test:
	@echo "=== PERFORMANCE TEST ==="
	@echo "Optimized compilation..."
	@$(MAKE) -s performance
	@echo "Quick test (modify params.txt for a longer test):"
	@time ./$(TARGET) 2>/dev/null | head -20 || echo "Test completed"