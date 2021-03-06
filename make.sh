#!/bin/sh

ROOT="$(dirname "$(readlink -f "$0")")"
EXENAME="gengeo"
CFLAGS="-std=c++14 -O2 -march=native -Wall -Wextra -Wno-unused-parameter -Wno-unused-function"

if [ x"$(uname -o)" = x"Msys" ]; then
	EXENAME="$EXENAME.exe"
fi

if [ "x$CXX" = "x" ]; then
	# Default to GCC
	CXX=g++
	# Prefer clang++ if available
	which clang++ > /dev/null 2> /dev/null
	if [ $? -eq 0 ]; then
		CXX=clang++
	fi
fi

set -x
$CXX $CFLAGS "$ROOT"/FastNoise/*.cpp "$ROOT"/src/*.cpp -o $EXENAME

