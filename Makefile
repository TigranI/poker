SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp, %.o, $(SRCS))

CXXFLAGS += -stdlib=libc++
CXXFLAGS += -std=c++11 -g  -I/usr/local/include

LDFLAGS += -L/usr/local/lib
LDFLAGS += -lnlopt_cxx -lm

CXX := clang++

EXEC := poker

all: $(EXEC)

$(EXEC) : $(OBJS)
		@echo Building recsystem, $(CXX) ...
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $@

%.o : %.cpp
		$(CXX) -c $(CXXFLAGS) $< -o $@
		@echo Building .o from .cpp $(CXX)

clean:
		rm $(EXEC) *.o

