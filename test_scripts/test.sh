#!/usr/bin/env bash

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <compress command> <source file>"
    exit 1
fi
GSC_COMMAND=../gsc
SOURCE_FILE="$1"
COMPRESS_FILE="$2"
DECOMPRESS_FILE="$3"

# if ! command -v "$COMPRESS_CMD" >/dev/null 2>&1; then
#     echo "Error: compress command '$COMPRESS_CMD' not found"
#     exit 1
# fi

function filesize_bytes {
    wc -c < "$1"
}

echo -e "\nCompressing file: $SOURCE_FILE"
original_size=$(filesize_bytes "$SOURCE_FILE")
start_time=$(date +%s.%N)

if [ "${SOURCE_FILE##*.}" == "bcf" ] ; then
    /usr/bin/time -f "%e %U %P %M" $GSC_COMMAND compress --bcf --all --out "$COMPRESS_FILE" "$SOURCE_FILE"
else
    /usr/bin/time -f "%e %U %P %M" $GSC_COMMAND compress --all --out "$COMPRESS_FILE" "$SOURCE_FILE"
fi

end_time=$(date +%s.%N)
compress_time=$(echo "scale=3; ($end_time - $start_time)" | bc)

compressed_size=$(echo "$(filesize_bytes "$COMPRESS_FILE.gti") + $(filesize_bytes "$COMPRESS_FILE.dbs")" | bc)
compression_ratio=$(echo "scale=3; $original_size / $compressed_size" | bc)
compression_speed=$(echo "scale=3; $original_size / ($compress_time * 1024 * 1024)" | bc)

echo -e "\nCompression result:"
echo "Original size: $original_size bytes"
echo "Compressed size: $compressed_size bytes"
printf "Compression ratio: %.3f\n" $compression_ratio
printf "Compression time: %.9f seconds\n" $compress_time
if [ "${SOURCE_FILE##*.}" == "gz" ] ; then
    printf "VCF.GZ file compression speed: %.2f MB/S\n" $compression_speed
elif [ "${SOURCE_FILE##*.}" == "bcf" ] ; then
    printf "BCF file compression speed: %.2f MB/S\n" $compression_speed
else
    printf "VCF file compression speed: %.2f MB/S\n" $compression_speed
fi


echo "***************************************************************"
echo -e "\nDecompressing file: $COMPRESS_FILE"
start_time=$(date +%s.%N)

/usr/bin/time -f "%e %U %P %M" $GSC_COMMAND decompress --bcf --all --out "$DECOMPRESS_FILE" "$COMPRESS_FILE"

end_time=$(date +%s.%N)
decompress_time=$(echo "scale=3; ($end_time - $start_time)" | bc)
original_size=$(filesize_bytes "$SOURCE_FILE")
decompress_size=$(filesize_bytes "$DECOMPRESS_FILE.bcf")
decompression_speed=$(echo "scale=3; $decompress_size / ($decompress_time * 1024 * 1024)" | bc)

echo -e "\nDecompression result:"
echo "Original size: $original_size bytes"
echo "BCF file Decompressed size: $decompress_size bytes"
printf "BCF file Decompression time: %.9f seconds\n" $decompress_time
printf "BCF file decompression speed: %.2f MB/S\n" $decompression_speed
echo "***************************************************************"

echo -e "\nChecking if compression/decompression works"
./check.sh "$SOURCE_FILE" "$DECOMPRESS_FILE.bcf"