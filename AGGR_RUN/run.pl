#!/usr/bin/perl
#
#
#
#
$pwd=`pwd`;
chomp $pwd;





#qsub -cwd -l  quick,h_vmem=200G
$cmd="qsub -cwd  -l quick,h_vmem=200G  -e $pwd/1.log -b y \"cellranger aggr --id=NEW2_1AGGR  --csv=new_1.agg.csv\" ";
print "$cmd\n";
system($cmd);


$cmd="qsub -cwd  -l quick,h_vmem=200G  -e $pwd/2.log -b y \"cellranger aggr --id=NEW2_2AGGR  --csv=new_2.agg.csv\" ";
print "$cmd\n";
system($cmd);



$cmd="qsub -cwd  -l quick,h_vmem=200G  -e $pwd/3.log -b y \"cellranger aggr --id=NEW2_3AGGR  --csv=new_3.agg.csv\" ";
print "$cmd\n";
system($cmd);






