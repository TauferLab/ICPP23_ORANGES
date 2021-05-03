#!/bin/bash

script_dir=$(pwd)
cd ..
find GRAAL -name "*.hpp" -o -name "*.cpp" -o -name "*.h" -o -name "*.c" | xargs clang-format -style=file -i
cd $script_dir
