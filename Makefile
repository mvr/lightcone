CXXC = clang++
CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3 -DNDEBUG -march=native -mtune=native -flto -fno-stack-protector -fomit-frame-pointer -fno-pic
# CXXFLAGS = -std=c++20 -Og -g3 -DDEBUG -Wall -Wextra -pedantic
LDFLAGS = -Wl,-stack_size -Wl,0x1000000

CXX = /opt/homebrew/opt/llvm/bin/clang++
# CXX = /usr/local/opt/llvm/bin/clang++
# LDFLAGS=-L/usr/local/opt/llvm/lib/c++ -Wl,-rpath,/usr/local/opt/llvm/lib/c++

CXXINCLUDE = -Itoml/include

all: LightCone

LightCone: LightCone.cpp *.hpp LifeAPI/*.hpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -o LightCone LightCone.cpp $(LDFLAGS)
