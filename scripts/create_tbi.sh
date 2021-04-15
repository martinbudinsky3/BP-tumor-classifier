#!/bin/bash

while getopts d:e: flag
do
	case "${flag}" in
		d) directory=${OPTARG};;
		e) end=${OPTARG};;			      
	esac
done
for d in "$directory"*/ ; do	
	for file in "$d"*"$end" ; do
		[ -f "$file" ] || break
		echo "$file"
		bgzip -c "$file" > "$file.gz"
		tabix -p vcf "$file.gz"
	done	
done
