#!/usr/bin/perl



$pwd=`pwd`;
chomp $pwd;



$cmd="qsub -cwd  -l quick,h_vmem=200G  -e $pwd/1.log -b y \"cellranger aggr --id=NEW2_1AGGR  --csv=new_1.agg.csv\" ";
print "$cmd\n";
system($cmd);
