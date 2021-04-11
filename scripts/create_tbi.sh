#!/bin/bash

for d in ../datasets/SpaceAndTime/VCFs/*/ ; do	
	for file in "$d"/*.vcf ; do
		echo "$file"
		#bgzip -c "$file" > "$file.gz"
		#tabix -p vcf "$file.gz"
	done	
done
