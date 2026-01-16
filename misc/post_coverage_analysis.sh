# Post-processing of Nextflow workflow output files

mkdir -p results

# Iterate over all *mosdepth.summary.txt files
for file in *.mosdepth.summary.txt; do
	# Extract sample ID from filename
	sid=${file%%.mosdepth.summary.txt}
	
	# Print the sample ID and the contents of the file
	# FS=OFS sets the field separator to tab
	awk -v id="${sid}" 'BEGIN { FS=OFS="\t" } { print $0, (NR==1 ? "sample_id" : id) }' "$file" > ${sid}_mosdepth.summary.txt
done

rm ./results/coverage.tsv
# Concatenate all *_mosdepth.summary.txt files
# The NR==1 condition removes the header line
awk 'FNR==1 && NR!=1 { next; } { print; }' *_mosdepth.summary.txt > ./results/coverage.tsv

# Replace "chrom" with "amplicon"
sed -i 's/chrom/amplicon/g' ./results/coverage.tsv

rm ./results/readcount.tsv
# Iterate over all *readcount.tsv files
for file in *readcount.tsv; do
	# Extract sample ID from filename
	sid=${file%%.readcount.tsv}
	
	# Print the sample ID and the contents of the file
	# FS=OFS sets the field separator to tab
	awk -v id="${sid}" 'BEGIN { FS=OFS="\t" } { print $0, id }' "$file" >> ./results/readcount.tsv
done

# Concatenate all *_reads_hnt.tsv files
# cat *reads_hnt.tsv > hnt_reads.tsv

# Replace "chrom" with "amplicon" and "reads" with "sample_id"
sed -i '1 i\amplicon	start	end	reads	sample_id' ./results/readcount.tsv

rm ./results/basecount.tsv

# Iterate over all *basecount.tsv files
for file in *basecount.tsv; do
	# Extract sample ID from filename
	sid=${file%%.basecount.tsv}
	
	# Print the sample ID and the contents of the file
	# FS=OFS sets the field separator to tab
	awk -v id="${sid}" 'BEGIN { FS=OFS="\t" } { print $0, id }' "$file" >> ./results/basecount.tsv
done

# Concatenate all *_reads_sm.tsv files
# cat *reads_sm.tsv > sm_reads.tsv

# replace "chrom" with "amplicon" and "reads" with "sample_id"
sed -i '1 i\amplicon	start	end	nbases	sample_id' ./results/basecount.tsv
