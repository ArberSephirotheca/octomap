# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++11 -Wall -Wextra -g

# Source files and executable name
SRCS := main.cpp OctreeBase.cpp
OBJS := $(SRCS:.cpp=.o)
EXEC := octomap

# Targets and rules
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)