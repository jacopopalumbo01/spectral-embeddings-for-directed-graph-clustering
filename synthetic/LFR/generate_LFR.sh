#!/bin/bash

# Get current directory
current_dir=${PWD##*/}               # to assign to a variable
current_dir=${current_dir:-/}        # to correct for the case where PWD is / (root)

right_dir="unweighted_directed"

if [[ "$current_dir" != "$right_dir" ]]; then
	echo "Current directory is \`$current_dir\`. The script must be run inside \`$right_dir\`"
	exit 1
fi

echo "Compiling the code..."
set -x
make -j
set +x

# Seeds to be used
seeds=(0 1 2 5 42 123 1234 1991 2001 2025)

# Mixing parameters
mu=(0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36 0.38 0.40)

# Create folder
mkdir graphs

for s in ${seeds[@]}; do
	
	echo "Generating LFR Benchmark Graph using seed=$s"
	echo "$s" >> time_seed.dat

	for m in ${mu[@]}; do

		echo "Using mixing parameter mu=$m"
		./benchmark -N 1000 -k 15 -maxk 50 -mu "$m" -minc 20 -maxc 50

		echo "Saving graph..."
		rm statistics.dat
		mv network.dat graphs/network-"$s"-"$m".dat
		mv community.dat graphs/community-"$s"-"$m".dat
	done
done
