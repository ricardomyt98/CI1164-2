#!/usr/bin/make

# folders
_SRC = ./src
_LIB = ./lib
_OBJ = ./obj

CC = gcc
RM = rm
CFLAGS = -lm -Wall -Werror -Wextra -I $(_LIB)

FILES =
FILES += partialDifferential
FILES += utils
FILES += pdeSolver

OBJECTS = $(addprefix $(_OBJ)/, $(addsuffix .o, $(FILES)))

$(_OBJ)/%.o: $(_SRC)/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

EXECUTABLE = pdeSolver

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) -o $(EXECUTABLE)

%(_OBJ)/%.o: $(_SRC)/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	$(RM) -f $(OBJECTS) $(EXECUTABLE)
