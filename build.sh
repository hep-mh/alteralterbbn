#! /usr/bin/env bash

clang src/*.c -ffast-math -std=c99 -Wall -Wpedantic -lm -O3 -o bin/alteralterbbn
