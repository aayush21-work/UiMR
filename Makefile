CXX      = g++
CXXFLAGS = -O3 -fopenmp -march=native
LDFLAGS  = -lm

SRC_DIR  = src
OBJ_DIR  = obj
TARGET   = benchmark

SRC      = $(SRC_DIR)/benchmark.cpp
OBJ      = $(OBJ_DIR)/benchmark.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

run: all
	./$(TARGET)

.PHONY: all clean run
