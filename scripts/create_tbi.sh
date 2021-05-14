#!/bin/bash

while getopts d:e:f: flag
do
	case "${flag}" in
		d) directory=${OPTARG};;
		e) end=${OPTARG};;
		f) in_file=${OPTARG};;		
	esac
done

create_tabix() {
	local file=$1
	bgzip -c "$file" > "$file.gz"
	tabix -p vcf "$file.gz"
}

if [ -z "$in_file" ]; 
	then
		for d in "$directory"/*/ ; do	
			for file in "$d"*"$end".vcf ; do
				[ -f "$file" ] || break
				echo "Creating tab index for $file"
				create_tabix $file
				echo "Tab index created"
			done	
		done;
	
	else
		echo "Creating tab index for  $in_file"
		create_tabix $in_file;
		echo "Tab index created"

fi


