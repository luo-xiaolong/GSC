#!/usr/bin/env bash
BCFTOOLS=./bcftools

cs1=$($BCFTOOLS view --no-version --threads 4 $1 | gawk ' !/^#/ {split( $8, a, ";" ); asort( a ); printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"; {for( i = 1; i < length(a); i++ ) printf("%s;", a[i]); printf(a[length(a)]"\t");} { for(i = 9; i<=NF; i++) printf($i"\t")} ;print ""}' | sed -e 's/0*e+/e+/g' | python3.8 xxh3.py)
cs2=$($BCFTOOLS view --no-version --threads 4 $2 | gawk ' !/^#/ {split( $8, a, ";" ); asort( a ); printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"; {for( i = 1; i < length(a); i++ ) printf("%s;", a[i]); printf(a[length(a)]"\t");} { for(i = 9; i<=NF; i++) printf($i"\t")} ;print ""}'  | python3.8 xxh3.py)

echo "original checksum:$cs1"
echo "restore checksum:$cs2"

if [ $cs1 == $cs2 ] ; then 
    echo OK; 
else 
    echo NOT OK; 
fi
