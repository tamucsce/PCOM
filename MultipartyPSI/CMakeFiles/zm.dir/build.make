# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ubuntu/MultipartyPSI

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ubuntu/MultipartyPSI

# Include any dependencies generated for this target.
include CMakeFiles/zm.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/zm.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/zm.dir/flags.make

CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.o: CMakeFiles/zm.dir/flags.make
CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.o: depends/ate-pairing/src/zm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/MultipartyPSI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.o -c /home/ubuntu/MultipartyPSI/depends/ate-pairing/src/zm.cpp

CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/MultipartyPSI/depends/ate-pairing/src/zm.cpp > CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.i

CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/MultipartyPSI/depends/ate-pairing/src/zm.cpp -o CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.s

CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.o: CMakeFiles/zm.dir/flags.make
CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.o: depends/ate-pairing/src/zm2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/MultipartyPSI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.o -c /home/ubuntu/MultipartyPSI/depends/ate-pairing/src/zm2.cpp

CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/MultipartyPSI/depends/ate-pairing/src/zm2.cpp > CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.i

CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/MultipartyPSI/depends/ate-pairing/src/zm2.cpp -o CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.s

# Object files for target zm
zm_OBJECTS = \
"CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.o" \
"CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.o"

# External object files for target zm
zm_EXTERNAL_OBJECTS =

libzm.a: CMakeFiles/zm.dir/depends/ate-pairing/src/zm.cpp.o
libzm.a: CMakeFiles/zm.dir/depends/ate-pairing/src/zm2.cpp.o
libzm.a: CMakeFiles/zm.dir/build.make
libzm.a: CMakeFiles/zm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/MultipartyPSI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libzm.a"
	$(CMAKE_COMMAND) -P CMakeFiles/zm.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/zm.dir/build: libzm.a

.PHONY : CMakeFiles/zm.dir/build

CMakeFiles/zm.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/zm.dir/cmake_clean.cmake
.PHONY : CMakeFiles/zm.dir/clean

CMakeFiles/zm.dir/depend:
	cd /home/ubuntu/MultipartyPSI && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/MultipartyPSI /home/ubuntu/MultipartyPSI /home/ubuntu/MultipartyPSI /home/ubuntu/MultipartyPSI /home/ubuntu/MultipartyPSI/CMakeFiles/zm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/zm.dir/depend

