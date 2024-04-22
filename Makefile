all: gsc

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not' 
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif

# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)

LIBS_DIR=lib
# LIBS_D= -L/usr/lib/x86_64-linux-gnu/ -lpython3.8 
INCLUDES_DIR=include 
# -I/usr/include/python3.8 -I/usr/local/lib/python3.8/dist-packages/numpy/core/include/
CFLAGS=-Wall -O3 -g -m64 -std=c++14 -pthread -I $(INCLUDES_DIR) -mpopcnt 
#CLINK=-O3 -lm -std=c++11 -lpthread -mpopcnt -lz
CLINK=-O3 -lm -std=c++14 -lpthread -mavx -mpopcnt -lz -lbz2 -llzma

ifeq ($(uname_S),Linux)
    CC=g++      
    # check if CPU supports SSE4.2 
    HAVE_SSE4=$(filter-out 0,$(shell grep sse4.2 /proc/cpuinfo | wc -l))
endif
ifeq ($(uname_S),Darwin)
    CC=clang++
    # check if CPU supports SSE4.2
    HAVE_SSE4=$(filter-out 0,$(shell  sysctl machdep.cpu.features| grep SSE4.2 - | wc -l))
endif

CFLAGS+=$(if $(HAVE_SSE4),-msse4.2)
ifeq ($(HAVE_SSE4),)
    CFLAGS+=-msse2
endif
-include src/*.d
.cpp.o:
	$(CC) $(CFLAGS) -c $< -MMD -o $@

gsc:	src/bit_memory.o \
	src/block_processing.o \
	src/vint_code.o	\
	src/decompression_reader.o \
	src/decompressor.o \
	src/compressor.o \
	src/main.o \
	src/samples.o \
	src/compression_reader.o \
	src/file_handle.o \
	src/utils.o \
	src/bsc.o \
	src/zstd_compress.o \
	include/cpp-mmf/memory_mapped_file.o
	$(CC) -o gsc \
	src/bit_memory.o \
	src/block_processing.o \
	src/vint_code.o	\
	src/decompression_reader.o \
	src/decompressor.o \
	src/compressor.o \
	src/main.o \
	src/samples.o \
	src/compression_reader.o \
	src/file_handle.o \
	src/utils.o \
	src/bsc.o \
	src/zstd_compress.o \
	include/cpp-mmf/memory_mapped_file.o \
	$(LIBS_DIR)/libhts.a \
	$(LIBS_DIR)/libsdsl.a \
	$(LIBS_DIR)/libbsc.a \
	$(LIBS_DIR)/libzstd.a \
	$(CLINK)


clean:
	-rm include/cpp-mmf/*.o
	-rm src/*.o
	-rm src/*.d
	-rm gsc
	
install:
	mkdir -p -m 755 $(exec_prefix)/bin
	cp gsc $(exec_prefix)/bin/
	
uninstall:
	rm  $(exec_prefix)/bin/gsc
	
