#!/usr/bin/env bash
mkdir xxhash_folder
cd xxhash_folder 
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh
./vcpkg integrate install
./vcpkg install xxhash