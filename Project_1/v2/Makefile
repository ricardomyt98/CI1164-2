# Folders
_DOC = doc
_SRC = src
_LIB = lib
_OBJ = obj
_LIK = likwidPerformance

# Programs
DOXYGEN = doxygen
CC = gcc
LINK = ln -rsf
RM = rm
FILE_RM = ${RM} -f
FOLDER_RM = ${RM} -rf

# Configuration variables
DOXYGEN_CONFIG = config
DOXYGEN_HTML = ${_DOC}/html
COMPILE_OBJ = -c
CFLAGS = -Wall -Ilib -std=gnu99 -lm $(OPTIMIZE_FLAGS) $(LIKWID_FLAGS)
LIKWID_FLAGS = -I/home/soft/likwid/include -L/home/soft/likwid/lib -I/usr/local/include -L/usr/local/lib -llikwid -DLIKWID_PERFMON
OPTIMIZE_FLAGS = -O3 -mavx -march=native
SRC_FILES = partialDifferential utils pdeSolver
OBJECTS = $(foreach src, $(SRC_FILES), ${_OBJ}/$(src).o)
EXEC = pdeSolver
DOXYGEN_COMPILED_FILES = ${DOXYGEN_HTML} ${_DOC}/latex
LIKWID_COMPILED_FILES = ${_LIK}/*


# Default rule
all: ${EXEC}

# Debug rule
debug: CFLAGS += -g -pg
debug: clean_files all

# Compilation of executable
${EXEC}: ${OBJECTS}
	${CC} ${OBJECTS} -o ${EXEC} ${CFLAGS}

# Doxygen documentation generation rule
doc: FORCE
	${DOXYGEN} ${_DOC}/${DOXYGEN_CONFIG}
	${LINK} ${DOXYGEN_HTML}/index.html Documentation.html

FORCE:

# Object compilation rule
${_OBJ}/%.o: ${_SRC}/%.c
	${CC} ${COMPILE_OBJ} $< -o $@ ${CFLAGS} ${DEBUG}

# Clean everything
clean: clean_files clean_doxygen clean_likwid

clean_files:
	${FILE_RM} ${OBJECTS} ${EXEC} gmon.out

clean_doxygen:
	${FOLDER_RM} ${DOXYGEN_COMPILED_FILES} Documentation.html

clean_likwid:
	${FOLDER_RM} ${LIKWID_COMPILED_FILES}
