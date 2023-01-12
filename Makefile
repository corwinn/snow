# Requirements: Open GL 2, FreeGLUT 3
MAKEFLAGS+=rR
CXX?=clang++
D=-fsanitize=undefined,leak,address -O0 -g
C=-std=c++14 -Wall -Wextra -Wshadow -Werror=shadow \
  -fno-rtti -fno-exceptions -fno-threadsafe-statics
L=-Wl,--as-needed -lGL -lglut
m:snow
%:%.cpp
	$(CXX) $(C) $< -o $@ $(L)

d:snow.cpp
	$(CXX) $(C) $(D) $< -o $@ $(L)

cake: universe