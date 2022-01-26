# MSMuTect-0.5-to-VCF
R code that can convert MSMuTect output file to .vcf

Usage:

$MSMuTect_FILE=my.MSMuTect.file

$INPUT_PAT=/my/input/folder

$OUTPUT_PATH=/my/output/folder

Rscript --vanilla MSMuTect_To_VCF.0.5.R $MSMuTect_FILE $INPUT_PATH $OUTPUT_PATH F 3 T T 5


7 argument must be supplied:

1-input file (char)

2-InputPath (char)

3-Output Path (char)

4-Single Sample Mode (Logical T/F)

5-Filtering mode (int 1..7):
                            1. Filter all non REF normals,
                            2. Filter normals homozygous to ALT,
                            3. Filter normals heterozygous to ALT,
                            4. No filters,
                            5. Exclude LOH: Remove variants which the normal is heterozygous and the tumor homozygous to one of the normal alleles.
                            6. Remove all normals with heterozygous samples.
                            7. Remove all normals with heterozygous or homozygous ALT alleles & all tumors that are not diploid heterozygous to REF.
                            
6-output-filtered MSMuTect file together with your vcf (Logical T/F)

7-Split multiAllelic loci in .vcf to single lines (Logical T/F)
