BUILD_DIR := build
CACHE := $(BUILD_DIR)/CMakeCache.txt

.PHONY: all test clean

all: build

build: $(CACHE)
	cmake --build $(BUILD_DIR)

$(CACHE):
	cmake -S . -B $(BUILD_DIR)

test: build
	ctest --test-dir $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

