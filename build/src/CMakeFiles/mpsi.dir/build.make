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
CMAKE_BINARY_DIR = /home/ubuntu/MultipartyPSI/build

# Include any dependencies generated for this target.
include src/CMakeFiles/mpsi.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/mpsi.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/mpsi.dir/flags.make

src/CMakeFiles/mpsi.dir/mpsi.cpp.o: src/CMakeFiles/mpsi.dir/flags.make
src/CMakeFiles/mpsi.dir/mpsi.cpp.o: ../src/mpsi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/MultipartyPSI/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/mpsi.dir/mpsi.cpp.o"
	cd /home/ubuntu/MultipartyPSI/build/src && /usr/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mpsi.dir/mpsi.cpp.o -c /home/ubuntu/MultipartyPSI/src/mpsi.cpp

src/CMakeFiles/mpsi.dir/mpsi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpsi.dir/mpsi.cpp.i"
	cd /home/ubuntu/MultipartyPSI/build/src && /usr/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/MultipartyPSI/src/mpsi.cpp > CMakeFiles/mpsi.dir/mpsi.cpp.i

src/CMakeFiles/mpsi.dir/mpsi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpsi.dir/mpsi.cpp.s"
	cd /home/ubuntu/MultipartyPSI/build/src && /usr/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/MultipartyPSI/src/mpsi.cpp -o CMakeFiles/mpsi.dir/mpsi.cpp.s

# Object files for target mpsi
mpsi_OBJECTS = \
"CMakeFiles/mpsi.dir/mpsi.cpp.o"

# External object files for target mpsi
mpsi_EXTERNAL_OBJECTS =

src/mpsi: src/CMakeFiles/mpsi.dir/mpsi.cpp.o
src/mpsi: src/CMakeFiles/mpsi.dir/build.make
src/mpsi: /usr/local/lib/libgmp.so
src/mpsi: /usr/local/lib/libgmpxx.so
src/mpsi: CMAKE_CURRENT_BINARY_DIR/libportable_binary_archive.a
src/mpsi: /usr/local/lib/libboost_program_options.a
src/mpsi: /usr/local/lib/libboost_system.a
src/mpsi: /usr/local/lib/libboost_log.a
src/mpsi: /usr/local/lib/libboost_serialization.a
src/mpsi: libzm.a
src/mpsi: /usr/local/lib/libgmp.so
src/mpsi: /usr/local/lib/libboost_chrono.a
src/mpsi: /usr/local/lib/libboost_filesystem.a
src/mpsi: /usr/local/lib/libboost_regex.a
src/mpsi: /usr/local/lib/libboost_thread.a
src/mpsi: /usr/local/lib/libboost_atomic.a
src/mpsi: src/CMakeFiles/mpsi.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/MultipartyPSI/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mpsi"
	cd /home/ubuntu/MultipartyPSI/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mpsi.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/mpsi.dir/build: src/mpsi

.PHONY : src/CMakeFiles/mpsi.dir/build

src/CMakeFiles/mpsi.dir/clean:
	cd /home/ubuntu/MultipartyPSI/build/src && $(CMAKE_COMMAND) -P CMakeFiles/mpsi.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/mpsi.dir/clean

src/CMakeFiles/mpsi.dir/depend:
	cd /home/ubuntu/MultipartyPSI/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/MultipartyPSI /home/ubuntu/MultipartyPSI/src /home/ubuntu/MultipartyPSI/build /home/ubuntu/MultipartyPSI/build/src /home/ubuntu/MultipartyPSI/build/src/CMakeFiles/mpsi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/mpsi.dir/depend

