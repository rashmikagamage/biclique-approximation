# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rashmika/biclique/sigmod_new/biclique-approximation/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rashmika/biclique/sigmod_new/biclique-approximation/src

# Include any dependencies generated for this target.
include CMakeFiles/exdensest.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/exdensest.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/exdensest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/exdensest.dir/flags.make

CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o: CMakeFiles/exdensest.dir/flags.make
CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o: runner/exactDensestSubgraph.cpp
CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o: CMakeFiles/exdensest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rashmika/biclique/sigmod_new/biclique-approximation/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o -MF CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o.d -o CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o -c /home/rashmika/biclique/sigmod_new/biclique-approximation/src/runner/exactDensestSubgraph.cpp

CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rashmika/biclique/sigmod_new/biclique-approximation/src/runner/exactDensestSubgraph.cpp > CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.i

CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rashmika/biclique/sigmod_new/biclique-approximation/src/runner/exactDensestSubgraph.cpp -o CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.s

# Object files for target exdensest
exdensest_OBJECTS = \
"CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o"

# External object files for target exdensest
exdensest_EXTERNAL_OBJECTS =

bin/exdensest: CMakeFiles/exdensest.dir/runner/exactDensestSubgraph.cpp.o
bin/exdensest: CMakeFiles/exdensest.dir/build.make
bin/exdensest: densestSubgraph/libexactFlowAlgorithm.a
bin/exdensest: CMakeFiles/exdensest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rashmika/biclique/sigmod_new/biclique-approximation/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/exdensest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/exdensest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/exdensest.dir/build: bin/exdensest
.PHONY : CMakeFiles/exdensest.dir/build

CMakeFiles/exdensest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/exdensest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/exdensest.dir/clean

CMakeFiles/exdensest.dir/depend:
	cd /home/rashmika/biclique/sigmod_new/biclique-approximation/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rashmika/biclique/sigmod_new/biclique-approximation/src /home/rashmika/biclique/sigmod_new/biclique-approximation/src /home/rashmika/biclique/sigmod_new/biclique-approximation/src /home/rashmika/biclique/sigmod_new/biclique-approximation/src /home/rashmika/biclique/sigmod_new/biclique-approximation/src/CMakeFiles/exdensest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/exdensest.dir/depend

