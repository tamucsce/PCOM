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
include CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/depend.make

# Include the progress variables for this target.
include CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/progress.make

# Include the compile flags for this target's objects.
include CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/flags.make

CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.o: CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/flags.make
CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.o: depends/LigeroLink/include/ligero/util/boost/portable_binary_iarchive.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/MultipartyPSI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.o"
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.o -c /home/ubuntu/MultipartyPSI/depends/LigeroLink/include/ligero/util/boost/portable_binary_iarchive.cpp

CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.i"
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/MultipartyPSI/depends/LigeroLink/include/ligero/util/boost/portable_binary_iarchive.cpp > CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.i

CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.s"
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/MultipartyPSI/depends/LigeroLink/include/ligero/util/boost/portable_binary_iarchive.cpp -o CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.s

CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.o: CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/flags.make
CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.o: depends/LigeroLink/include/ligero/util/boost/portable_binary_oarchive.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/MultipartyPSI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.o"
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.o -c /home/ubuntu/MultipartyPSI/depends/LigeroLink/include/ligero/util/boost/portable_binary_oarchive.cpp

CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.i"
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/MultipartyPSI/depends/LigeroLink/include/ligero/util/boost/portable_binary_oarchive.cpp > CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.i

CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.s"
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/MultipartyPSI/depends/LigeroLink/include/ligero/util/boost/portable_binary_oarchive.cpp -o CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.s

# Object files for target portable_binary_archive
portable_binary_archive_OBJECTS = \
"CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.o" \
"CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.o"

# External object files for target portable_binary_archive
portable_binary_archive_EXTERNAL_OBJECTS =

CMAKE_CURRENT_BINARY_DIR/libportable_binary_archive.a: CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_iarchive.cpp.o
CMAKE_CURRENT_BINARY_DIR/libportable_binary_archive.a: CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/include/ligero/util/boost/portable_binary_oarchive.cpp.o
CMAKE_CURRENT_BINARY_DIR/libportable_binary_archive.a: CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/build.make
CMAKE_CURRENT_BINARY_DIR/libportable_binary_archive.a: CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/MultipartyPSI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libportable_binary_archive.a"
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && $(CMAKE_COMMAND) -P CMakeFiles/portable_binary_archive.dir/cmake_clean_target.cmake
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/portable_binary_archive.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/build: CMAKE_CURRENT_BINARY_DIR/libportable_binary_archive.a

.PHONY : CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/build

CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/clean:
	cd /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR && $(CMAKE_COMMAND) -P CMakeFiles/portable_binary_archive.dir/cmake_clean.cmake
.PHONY : CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/clean

CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/depend:
	cd /home/ubuntu/MultipartyPSI && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/MultipartyPSI /home/ubuntu/MultipartyPSI/depends/LigeroLink /home/ubuntu/MultipartyPSI /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR /home/ubuntu/MultipartyPSI/CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMAKE_CURRENT_BINARY_DIR/CMakeFiles/portable_binary_archive.dir/depend
