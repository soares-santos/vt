# CMake generated Testfile for 
# Source directory: /home/fsobreira/athena/athena_1.7
# Build directory: /home/fsobreira/athena/athena_1.7/build
# 
# This file replicates the SUBDIRS() and ADD_TEST() commands from the source
# tree CMakeLists.txt file, skipping any SUBDIRS() or ADD_TEST() commands
# that are excluded by CMake control structures, i.e. IF() commands.
ADD_TEST(test_suite "/home/fsobreira/athena/athena_1.7/bin/test_suite_athena.py" "-t" "15" "-P" "/home/fsobreira/athena/athena_1.7/bin" "-T" "/home/fsobreira/athena/athena_1.7/test")
SUBDIRS(src)
