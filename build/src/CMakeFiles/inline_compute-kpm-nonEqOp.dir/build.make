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
CMAKE_COMMAND = /home/santiago/miniconda3/bin/cmake

# The command to remove a file.
RM = /home/santiago/miniconda3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build

# Include any dependencies generated for this target.
include src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/flags.make

src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o: src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/flags.make
src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o: /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/src/inline_compute-kpm-nonEqOp.cpp
src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o: src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o"
	cd /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/src && /opt/intel/oneapi/compiler/2024.0/bin/icpx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o -MF CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o.d -o CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o -c /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/src/inline_compute-kpm-nonEqOp.cpp

src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.i"
	cd /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/src && /opt/intel/oneapi/compiler/2024.0/bin/icpx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/src/inline_compute-kpm-nonEqOp.cpp > CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.i

src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.s"
	cd /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/src && /opt/intel/oneapi/compiler/2024.0/bin/icpx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/src/inline_compute-kpm-nonEqOp.cpp -o CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.s

# Object files for target inline_compute-kpm-nonEqOp
inline_compute__kpm__nonEqOp_OBJECTS = \
"CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o"

# External object files for target inline_compute-kpm-nonEqOp
inline_compute__kpm__nonEqOp_EXTERNAL_OBJECTS =

inline_compute-kpm-nonEqOp: src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/inline_compute-kpm-nonEqOp.o
inline_compute-kpm-nonEqOp: src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/build.make
inline_compute-kpm-nonEqOp: lib/libkpm_lib.a
inline_compute-kpm-nonEqOp: src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../inline_compute-kpm-nonEqOp"
	cd /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/inline_compute-kpm-nonEqOp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/build: inline_compute-kpm-nonEqOp
.PHONY : src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/build

src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/clean:
	cd /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/src && $(CMAKE_COMMAND) -P CMakeFiles/inline_compute-kpm-nonEqOp.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/clean

src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/depend:
	cd /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/src /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/src /home/santiago/Documents/ICN2/Codes/linqt-2.0.0_beta/build/src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/inline_compute-kpm-nonEqOp.dir/depend

