# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.24.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.24.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/nicolas/Desktop/frechet_distance

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/nicolas/Desktop/frechet_distance/build

# Include any dependencies generated for this target.
include CMakeFiles/create_benchmark_decider.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/create_benchmark_decider.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/create_benchmark_decider.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/create_benchmark_decider.dir/flags.make

CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o: CMakeFiles/create_benchmark_decider.dir/flags.make
CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/create_benchmark_decider.cpp
CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o: CMakeFiles/create_benchmark_decider.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o -MF CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o.d -o CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/create_benchmark_decider.cpp

CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/create_benchmark_decider.cpp > CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.i

CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/create_benchmark_decider.cpp -o CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.s

# Object files for target create_benchmark_decider
create_benchmark_decider_OBJECTS = \
"CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o"

# External object files for target create_benchmark_decider
create_benchmark_decider_EXTERNAL_OBJECTS = \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/frechet_light.cpp.o" \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/frechet_naive.cpp.o" \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/geometry_basics.cpp.o" \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/filter.cpp.o" \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/orth_range_search.cpp.o" \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/parser.cpp.o" \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/query.cpp.o" \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/times.cpp.o" \
"/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/common.dir/src/curve.cpp.o"

create_benchmark_decider: CMakeFiles/create_benchmark_decider.dir/src/create_benchmark_decider.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/frechet_light.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/frechet_naive.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/geometry_basics.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/filter.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/orth_range_search.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/parser.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/query.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/times.cpp.o
create_benchmark_decider: CMakeFiles/common.dir/src/curve.cpp.o
create_benchmark_decider: CMakeFiles/create_benchmark_decider.dir/build.make
create_benchmark_decider: CMakeFiles/create_benchmark_decider.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable create_benchmark_decider"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/create_benchmark_decider.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/create_benchmark_decider.dir/build: create_benchmark_decider
.PHONY : CMakeFiles/create_benchmark_decider.dir/build

CMakeFiles/create_benchmark_decider.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/create_benchmark_decider.dir/cmake_clean.cmake
.PHONY : CMakeFiles/create_benchmark_decider.dir/clean

CMakeFiles/create_benchmark_decider.dir/depend:
	cd /Users/nicolas/Desktop/frechet_distance/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nicolas/Desktop/frechet_distance /Users/nicolas/Desktop/frechet_distance /Users/nicolas/Desktop/frechet_distance/build /Users/nicolas/Desktop/frechet_distance/build /Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/create_benchmark_decider.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/create_benchmark_decider.dir/depend

