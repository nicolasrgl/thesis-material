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
include CMakeFiles/export_freespace_diagram.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/export_freespace_diagram.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/export_freespace_diagram.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/export_freespace_diagram.dir/flags.make

CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/export_freespace_diagram.cpp
CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/export_freespace_diagram.cpp

CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/export_freespace_diagram.cpp > CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/export_freespace_diagram.cpp -o CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/frechet_light.cpp
CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/frechet_light.cpp

CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/frechet_light.cpp > CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/frechet_light.cpp -o CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/frechet_naive.cpp
CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/frechet_naive.cpp

CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/frechet_naive.cpp > CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/frechet_naive.cpp -o CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/geometry_basics.cpp
CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/geometry_basics.cpp

CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/geometry_basics.cpp > CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/geometry_basics.cpp -o CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/filter.cpp
CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/filter.cpp

CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/filter.cpp > CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/filter.cpp -o CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/freespace_light_vis.cpp
CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/freespace_light_vis.cpp

CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/freespace_light_vis.cpp > CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/freespace_light_vis.cpp -o CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/orth_range_search.cpp
CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/orth_range_search.cpp

CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/orth_range_search.cpp > CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/orth_range_search.cpp -o CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/parser.cpp
CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/parser.cpp

CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/parser.cpp > CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/parser.cpp -o CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/query.cpp
CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/query.cpp

CMakeFiles/export_freespace_diagram.dir/src/query.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/query.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/query.cpp > CMakeFiles/export_freespace_diagram.dir/src/query.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/query.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/query.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/query.cpp -o CMakeFiles/export_freespace_diagram.dir/src/query.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/times.cpp
CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/times.cpp

CMakeFiles/export_freespace_diagram.dir/src/times.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/times.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/times.cpp > CMakeFiles/export_freespace_diagram.dir/src/times.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/times.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/times.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/times.cpp -o CMakeFiles/export_freespace_diagram.dir/src/times.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/certificate.cpp
CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/certificate.cpp

CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/certificate.cpp > CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/certificate.cpp -o CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.s

CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o: CMakeFiles/export_freespace_diagram.dir/flags.make
CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o: /Users/nicolas/Desktop/frechet_distance/src/curve.cpp
CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o: CMakeFiles/export_freespace_diagram.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o -MF CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o.d -o CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o -c /Users/nicolas/Desktop/frechet_distance/src/curve.cpp

CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nicolas/Desktop/frechet_distance/src/curve.cpp > CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.i

CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nicolas/Desktop/frechet_distance/src/curve.cpp -o CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.s

# Object files for target export_freespace_diagram
export_freespace_diagram_OBJECTS = \
"CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o" \
"CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o"

# External object files for target export_freespace_diagram
export_freespace_diagram_EXTERNAL_OBJECTS =

export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/export_freespace_diagram.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/frechet_light.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/frechet_naive.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/geometry_basics.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/filter.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/freespace_light_vis.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/orth_range_search.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/parser.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/query.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/times.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/certificate.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/src/curve.cpp.o
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/build.make
export_freespace_diagram: CMakeFiles/export_freespace_diagram.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/nicolas/Desktop/frechet_distance/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX executable export_freespace_diagram"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/export_freespace_diagram.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/export_freespace_diagram.dir/build: export_freespace_diagram
.PHONY : CMakeFiles/export_freespace_diagram.dir/build

CMakeFiles/export_freespace_diagram.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/export_freespace_diagram.dir/cmake_clean.cmake
.PHONY : CMakeFiles/export_freespace_diagram.dir/clean

CMakeFiles/export_freespace_diagram.dir/depend:
	cd /Users/nicolas/Desktop/frechet_distance/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nicolas/Desktop/frechet_distance /Users/nicolas/Desktop/frechet_distance /Users/nicolas/Desktop/frechet_distance/build /Users/nicolas/Desktop/frechet_distance/build /Users/nicolas/Desktop/frechet_distance/build/CMakeFiles/export_freespace_diagram.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/export_freespace_diagram.dir/depend

