#!/usr/bin/perl




 open my $fh, '<', 'region_genotypeMatrix.traw' or die "Cannot open file: $!";

 # Read and print the header
 my $header = <$fh>;
# print $header;
print "CHR\tPOS\n";
 # Iterate through the lines in the file
 while (my $line = <$fh>) {
     chomp $line;
next if $line =~ m/NA/g;
next if $line =~ m/KI27/g;
         # Split the line into columns
             my @columns = split /\t/, $line;

                 # Extract the values of 0_T124 and 0_T125 columns
                     my $value_0_T124 = $columns[6];
                         my $value_0_T125 = $columns[7];
          #  print "COL:$columns[7]\t$columns[8]\n";
                             # Print the line if there is a difference in 0_T124 and 0_T125
                            next if $columns[6] == 0;
			    next if $columns[7] == 0; 
                                 if ($value_0_T124 != $value_0_T125) {
                                        print $line, "\n";
                                        #print "$columns[0]\t$columns[3]\t$columns[3]\n"; 
                                       #chrCHR	POS
                                       #
                                       #print "$columns[0]\t$columns[3]\n"; 
                                             }
                                             }

                                             # Close the file
                                             close $fh;
