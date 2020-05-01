#!/usr/bin/perl
use warnings;

#This program extracts the flanking reads around the insertion break point and also their mates mapping on mitochondrial DNA from all the samples the event was discovered in. Then dynamically creates output files with the reads.  
#Usage : BrkPT_read_extract.pl <File with Numts and samples they are found in bedFormat.txt> <directory with BAM files>

$number  = $#ARGV;
if ($#ARGV < 1)
{
	print "please enter it in the following format: $0<input file name><directory_with_BAMfiles>\n";
	exit 1;
}

$events = $ARGV[0];
$dir = $ARGV[1];



open( fname1, $events )     or die("error opening file $events\n");


while ( $line1 = <fname1> ) 
{
    chomp($line1);
    @val_array = split( /\t/, $line1 );

	my $chr = $val_array[0];
    my $brk_pt1 = $val_array[1];
    my $brk_pt2 = $val_array[2];
    my $sample = $val_array[3];
    
	my $brk_pt3 = $brk_pt1 + 500;
    my $brk_pt4 = $brk_pt2 - 500;
    
	open (foutname1, ">$chr\_$brk_pt1.del") or die ("error opening file $sam\n");
    
    my @sam_val = split ( ',', $sample );
    
    my $sam_value = @sam_val;

    for ( my $i = 0; $i < $sam_value; $i++) 
    	{
    
        	my $sam = $sam_val[$i];
        
        
    		open( SAM, "samtools view $dir$sam.bam |" ) || die "Could not open $sam for input, $!\n";
       
			while (<SAM>) 
			{
    			chomp;
    			my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, $opt ) = split(/\t/);
 
        		if ($pos >= $brk_pt1 && $pos <= $brk_pt2 && $rname eq $chr)
        		{
            		if ($rnext eq 'chrM')
            		{
                		print  foutname1 ">$rname\t$pos\t$rnext\t$pnext\n$seq\n";
            		}
            		elsif ($pos >= $brk_pt3 && $pos <= $brk_pt4 && $rnext eq '=' )
            		{
                		if ($cigar =~ m/S/)
                		{
                    		print foutname1 ">$rname\t$pos\t$rnext\t$pnext\t$cigar\n$seq\n";
                		}
            		}
        		}
        
        		elsif ($rname eq 'chrM' )
        		{
            		if ($rnext eq $chr && $pnext >= $brk_pt1 && $pnext <= $brk_pt2)
            		{
                		print foutname1 ">$rname\t$pos\t$rnext\t$pnext\n$seq\n";
            		}
        		}
        
    		}

	close (SAM);

		}
close (foutname1);
}
