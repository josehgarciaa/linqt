# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_SOURCE_DIR = /home/jgarcia/linqt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jgarcia/linqt/build

# Include any dependencies generated for this target.
include src/CMakeFiles/kpm_lib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/kpm_lib.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/kpm_lib.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/kpm_lib.dir/flags.make

src/CMakeFiles/kpm_lib.dir/dummy.cpp.o: src/CMakeFiles/kpm_lib.dir/flags.make
src/CMakeFiles/kpm_lib.dir/dummy.cpp.o: ../src/dummy.cpp
src/CMakeFiles/kpm_lib.dir/dummy.cpp.o: src/CMakeFiles/kpm_lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jgarcia/linqt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/kpm_lib.dir/dummy.cpp.o"
	cd /home/jgarcia/linqt/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/kpm_lib.dir/dummy.cpp.o -MF CMakeFiles/kpm_lib.dir/dummy.cpp.o.d -o CMakeFiles/kpm_lib.dir/dummy.cpp.o -c /home/jgarcia/linqt/src/dummy.cpp

src/CMakeFiles/kpm_lib.dir/dummy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/kpm_lib.dir/dummy.cpp.i"
	cd /home/jgarcia/linqt/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jgarcia/linqt/src/dummy.cpp > CMakeFiles/kpm_lib.dir/dummy.cpp.i

src/CMakeFiles/kpm_lib.dir/dummy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/kpm_lib.dir/dummy.cpp.s"
	cd /home/jgarcia/linqt/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jgarcia/linqt/src/dummy.cpp -o CMakeFiles/kpm_lib.dir/dummy.cpp.s

src/CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o: src/CMakeFiles/kpm_lib.dir/flags.make
src/CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o: ../src/chebyshev_moments1D.cpp
src/CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o: src/CMakeFiles/kpm_lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jgarcia/linqt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o"
	cd /home/jgarcia/linqt/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o -MF CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o.d -o CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o -c /home/jgarcia/linqt/src/chebyshev_moments1D.cpp

src/CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.i"
	cd /home/jgarcia/linqt/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jgarcia/linqt/src/chebyshev_moments1D.cpp > CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.i

src/CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.s"
	cd /home/jgarcia/linqt/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jgarcia/linqt/src/chebyshev_moments1D.cpp -o CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.s

# Object files for target kpm_lib
kpm_lib_OBJECTS = \
"CMakeFiles/kpm_lib.dir/dummy.cpp.o" \
"CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o"

# External object files for target kpm_lib
kpm_lib_EXTERNAL_OBJECTS =

lib/libkpm_lib.a: src/CMakeFiles/kpm_lib.dir/dummy.cpp.o
lib/libkpm_lib.a: src/CMakeFiles/kpm_lib.dir/chebyshev_moments1D.cpp.o
lib/libkpm_lib.a: src/CMakeFiles/kpm_lib.dir/build.make
lib/libkpm_lib.a: src/CMakeFiles/kpm_lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jgarcia/linqt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library ../lib/libkpm_lib.a"
	cd /home/jgarcia/linqt/build/src && $(CMAKE_COMMAND) -P CMakeFiles/kpm_lib.dir/cmake_clean_target.cmake
	cd /home/jgarcia/linqt/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/kpm_lib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/kpm_lib.dir/build: lib/libkpm_lib.a
.PHONY : src/CMakeFiles/kpm_lib.dir/build

src/CMakeFiles/kpm_lib.dir/clean:
	cd /home/jgarcia/linqt/build/src && $(CMAKE_COMMAND) -P CMakeFiles/kpm_lib.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/kpm_lib.dir/clean

src/CMakeFiles/kpm_lib.dir/depend:
	cd /home/jgarcia/linqt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jgarcia/linqt /home/jgarcia/linqt/src /home/jgarcia/linqt/build /home/jgarcia/linqt/build/src /home/jgarcia/linqt/build/src/CMakeFiles/kpm_lib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/kpm_lib.dir/depend

