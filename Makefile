CXX = clang++
CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3 -DNDEBUG -march=native -mtune=native -flto -fno-stack-protector -fomit-frame-pointer -fno-pic
# CXXFLAGS = -std=c++20 -Og -g3 -DDEBUG -Wall -Wextra -pedantic
LDFLAGS =

# CXX = /opt/homebrew/opt/llvm/bin/clang++
# CXX = /usr/local/opt/llvm/bin/clang++
# LDFLAGS=-L/usr/local/opt/llvm/lib/c++ -Wl,-rpath,/usr/local/opt/llvm/lib/c++

CXXINCLUDE = -Itoml/include

all: Lightcone

Lightcone: Lightcone.cpp *.hpp LifeAPI/*.hpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -o Lightcone Lightcone.cpp $(LDFLAGS)
