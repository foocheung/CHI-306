#!/usr/bin/perl

# Command 1: Convert VCF to binary PLINK format
system("plink2 --vcf twins_liftover.vcf --make-bed --out region --allow-extra-chr --max-alleles 2");

# Command 2: Recode and transpose PLINK binary files
system("plink2 --bfile region --recode A-transpose --out region_genotypeMatrix --allow-extra-chr");

# Command 3: Parse TRAW file (assuming parse_traw.pl is a Perl script)
system("parse_traw.pl > filter");

# Command 4: Filter VCF using positions from the parsed TRAW file
system("vcftools --vcf twins_liftover.vcf --positions filter --recode --out V2_twins_liftover.vcf");

# Command 5: Compress the filtered VCF using bgzip
system("bgzip V2_twins_liftover.vcf");

# Command 6: Index the compressed VCF using tabix
system("tabix -f V2_twins_liftover.vcf.gz");

# Command 7: Index the original VCF (assuming nontwins_liftover.vcf.gz is the original VCF)
system("tabix -f nontwins_liftover.vcf.gz");

# Command 8: Merge the indexed VCF files using bcftools
system("bcftools merge nontwins_liftover.vcf.gz V2_twins_liftover.vcf.recode.vcf.gz > new.vcf");



#The script assumes that the required tools (plink2, parse_traw.pl, vcftools, tabix, and bcftools) are in your system's PATH
