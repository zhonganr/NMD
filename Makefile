CC		:= g++
CFLAGS	:= -std=c++17 -Wall -Wextra -g

BIN		:= bin
SRC		:= src
INCLUDE	:= include "D:\LLVM\include" "D:\LLVM\include\QtWidgets" "D:\LLVM\include\QtGui" "D:\LLVM\include\QtCore"
LIB		:= lib "D:\LLVM\lib" 

LIBRARIES	:= glfw3.dll -lfreetype -lz -lQt5Core -lQt5Gui -lQt5Widgets

ifeq ($(OS),Windows_NT)
EXECUTABLE	:= Magnetoelastic_Tool.exe
SOURCEDIRS	:= $(SRC)
INCLUDEDIRS	:= $(INCLUDE)
LIBDIRS		:= $(LIB)
else
EXECUTABLE	:= Magnetoelastic_Tool
SOURCEDIRS	:= $(shell find $(SRC) -type d)
INCLUDEDIRS	:= $(shell find $(INCLUDE) -type d)
LIBDIRS		:= $(shell find $(LIB) -type d)
endif

CINCLUDES	:= $(patsubst %,-I%, $(INCLUDEDIRS:%/=%))
CLIBS		:= $(patsubst %,-L%, $(LIBDIRS:%/=%))

SOURCES_CXX	:= $(wildcard $(patsubst %,%/*.cpp, $(SOURCEDIRS)))
SOURCES_C	:= $(wildcard $(patsubst %,%/*.c, $(SOURCEDIRS)))
#OBJECTS		:= $(SOURCES:.cpp=.o)
OBJECTS		:= $(SOURCES_CXX) $(SOURCES_C)

all: $(BIN)/$(EXECUTABLE)

.PHONY: clean
clean:
	-$(RM) $(BIN)/$(EXECUTABLE)
#	-$(RM) $(OBJECTS)

run: all
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(CINCLUDES) $(CLIBS) $(LIBRARIES) -o $@