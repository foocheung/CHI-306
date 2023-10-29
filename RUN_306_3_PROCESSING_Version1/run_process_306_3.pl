#!/usr/bin/perl


$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=ALL --lanes=1,2,3,4,5,6,7,8";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane1245_CD19+S1-  --lanes=1,2,4,5";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane3_S1+  --lanes=3";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane7_CD95+CD4+  --lanes=7";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane8_CD95+CD8+  --lanes=8";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane6_S1-  --lanes=6";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane1_CD19+S1-  --lanes=1";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane2_CD19+S1-  --lanes=2";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane4_CD19+S1-  --lanes=4";
system($cmd);
print "$cmd\n";

$cmd="Rscript ./process_lanes_v03.R --rawdir=./ --demuxdir=./DEMUX/ --save_name=lane5_CD19+S1-  --lanes=5";
system($cmd);
print "$cmd\n";
