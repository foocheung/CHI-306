#!/usr/bin/perl

$pwd=`pwd`;
chomp $pwd;


system("module load seqkit/2.3.1");

system("mkdir DOWNSAMPLE_BY_READS");
while(<>){
chomp;
next if $. ==1;

($a,$b,$count,$dir)=split(/\t/,$_);

$gex=$a;
$gex=~ s/^.*\_//g;

system("mkdir DOWNSAMPLE_BY_READS/GEX$gex");

$cmd="qsub -cwd -l  quick,h_vmem=200G -e $pwd/DOWNSAMPLE_BY_READS/$gex\.log -o $pwd/DOWNSAMPLE_BY_READS/$gex\.qsub -b y \"seqkit sample --two-pass -s100 -n $b GEX$gex/23_306_3_GEX$gex\_S$count\_R1_001.fastq.gz  -o DOWNSAMPLE_BY_READS/GEX$gex/N
EW_S23_306_3_GEX$gex\_S$count\_R1_001.fastq.gz\" ";
print "$cmd\n";
system($cmd);



$cmd="qsub -cwd -l  quick,h_vmem=200G -e $pwd/DOWNSAMPLE_BY_READS/$gex\.log -o $pwd/DOWNSAMPLE_BY_READS/$gex\.qsub -b y \"seqkit sample --two-pass -s100 -n  $b GEX$gex/23_306_3_GEX$gex\_S$count\_R2_001.fastq.gz  -o DOWNSAMPLE_BY_READS/GEX$gex/
NEW_S23_306_3_GEX$gex\_S$count\_R2_001.fastq.gz\" ";
print "$cmd\n";
system($cmd);

}
