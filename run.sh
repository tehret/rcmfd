#!/usr/bin/env bash
set -e

# Parameters
gs=$1
ps=$2
auto=$3
tau=$4
eps=$5

# Do the forgery detection
if [ "$auto" = True ]; then
    main -im input_0.png -ps $ps -gs $gs -auto true -tau $eps -vo output.png
else
    main -im input_0.png -ps $ps -gs $gs -auto false -tau $tau -vo output.png
fi;

# Check if something has been detected
matches=`wc -l < data_matches.csv`

# Produce the final results
if [ "$matches" -gt 1 ]; then
    echo "isforgery=1" > algo_info.txt
    echo "This image IS a forgery" > result.txt
else #[ "$Type" == 2 ]; then
    echo "This image IS NOT a forgery" > result.txt
fi;
