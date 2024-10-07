# Compiler and flags
CXX = g++
CXXFLAGS = -I /home/duc/miniconda3/envs/py310/include/python3.10
LDFLAGS = -lyaml-cpp -lpython3.10 -lgfortran
TARGET = KamodoWrapperExample

USER_FORTRAN_FILE = main.f90
OBJECTS = KamodoWrapper.o FortranWrapper.o main.o

# Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

KamodoWrapper.o: KamodoWrapper.cpp
	$(CXX) $(CXXFLAGS) $< -c $@

FortranWrapper.o: FortranWrapper.f90
	$(CXX) $< -c $@

main.o: $(USER_FORTRAN_FILE) FortranWrapper.o KamodoWrapper.o
	$(CXX) $(LDFLAGS) $< -c $@

# Clean rule
clean:
	rm -f $(OBJECTS) $(TARGET) *.mod

.PHONY: all clean
