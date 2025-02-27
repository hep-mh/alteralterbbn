# Specify the compiler
CC ?= clang

# Specify the compiler options
OPTS = -ffast-math -std=c99 -Wall -Wpedantic -lm -O3

default: alteralterbbn

alteralterbbn: src/*.c
	$(CC) $^ $(OPTS) -o bin/$@