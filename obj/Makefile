CC := g++-7
CCFLAGS := -std=c++14
INCLUDE := ../src
ROOT_FLAGS := $(shell echo $$(root-config --cflags --libs) )
ROOT_FLAGS := $(shell echo $(ROOT_FLAGS) | sed 's/-std=c++11//;s/-stdlib=libc++//g' )

VPATH = ../src:../main
SOURCES := $(wildcard ../src/*.cpp ../main/*.cpp)
OBJECTS := $(notdir $(SOURCES:.cpp=.o))

# OBJECTS := $(addprefix ../bin/, $(notdir $(SOURCES:.cpp=.o)))


#####################################################

all: ../bin/Projecto.exe

../bin/Projecto.exe: $(OBJECTS)
	@echo $(OBJECTS)
	$(CC) -o $@ $^ -I $(INCLUDE) $(ROOT_FLAGS)
	# ./bin/Projecto.exe

./%.o: %.cpp
	$(CC) $(CCFLAGS) -c $^ -o $@ -I $(INCLUDE) $(ROOT_FLAGS)

clean:
	rm ./*.o ../bin/Projecto.exe
