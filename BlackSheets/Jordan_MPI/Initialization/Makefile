CXX = mpic++
CXXFLAGS = -isystem /usr/lib/x86_64-linux-gnu/openmpi/include \
	-O3 -mfpmath=sse -fstack-protector-all -g \
	-Wall -Wextra -Werror -pedantic -pedantic-errors \
        -Wformat-security -Wmissing-format-attribute \
        -Wwrite-strings -Wcast-align -Wno-long-long \
        -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual

SRCS = file.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = program

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp func.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
