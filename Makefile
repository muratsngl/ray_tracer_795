# Define the C++ compiler to use
CXX = g++

# Define directories
SRCDIR = ./source
HDRDIR = ./header
OBJDIR = ./obj
THIRDPARTY = ./third_party

# Define compiler flags
# -std=c++17: Use the C++17 standard
# -Wall -Wextra: Enable most warnings
# -g: Include debug symbols
# -O3: Maximum optimization level
# -march=native: Enable all instruction subsets supported by the local machine
# -mavx2: Enable AVX2 instructions (comment out if your CPU doesn't support it)
# -mfma: Enable FMA (Fused Multiply-Add) instructions
# -fopenmp: Enable OpenMP support for parallel processing
# -I$(HDRDIR): Tell the compiler to look for headers in the ./header directory
# -I$(THIRDPARTY)/include: Tell the compiler to look for third-party headers (xsimd)
CXXFLAGS = -std=c++17 -Wall -Wextra -g -O2 -march=native -mavx2 -mfma -fopenmp -I$(HDRDIR) -I$(THIRDPARTY)/include

# Define the name of the final executable
TARGET = raytracer

# --- Automatic File Discovery ---
# Find all .cpp files in the source directory
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
# Create a list of object files that will be placed in the obj directory
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
# Find all header files that should trigger recompilation
HEADERS = $(wildcard $(HDRDIR)/*.hpp)

# --- Build Rules ---

# The default rule, executed when you type 'make'
.PHONY: all
all: $(TARGET)

# Rule to link the final executable from the object files
$(TARGET): $(OBJECTS)
	@echo "==> Linking..."
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "==> Build finished: ./"$(TARGET)

# Rule to compile a source file (.cpp) into an object file (.o)
# This rule creates the object directory if it doesn't exist first
# Now depends on ALL header files, so any header change triggers recompilation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	@mkdir -p $(OBJDIR)
	@echo "==> Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to clean up the project (remove executable and object files)
.PHONY: clean
clean:
	@echo "==> Cleaning up build files..."
	rm -rf $(OBJDIR) $(TARGET)