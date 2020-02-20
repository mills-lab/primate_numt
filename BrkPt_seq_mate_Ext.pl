#!/usr/bin/perl


my $dir = "/mnt/EXT/Mills-scratch/reference/primate/orang/";
my $events = "/mnt/EXT/Mills-data/gdayama/numts/primates/SRA_dinumt/analysis/brkPT/orang_Numt_final.txt";
my $fout1 = "/mnt/EXT/Mills-data/gdayama/numts/primates/SRA_dinumt/analysis/brkPT/orang_Numt_flanks.cmds";


open( fname1, $events )     or die("error opening file $events\n");



while ( $line1 = <fname1> ) {
    chomp($line1);
    @val_array = split( /\t/, $line1 );

	my $chr = $val_array[0];
    my $brk_pt1 = $val_array[1];
    my $brk_pt2 = $val_array[2];
    my $sample = $val_array[3];
    
    while ( $sample ) {
    chomp($sample);
    @sam_val = split ( /,/, $sample ); 
    
    open (foutname1, ">$sam_val") or die ("error opening file $fout1\n");
    open( SAM, "samtools view $sam_val* |" ) || die "Could not open $sam_val for input, $!\n";
	while (<SAM>) {
    chomp;
    my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, $opt ) = split(/\t/);
		
		if ($rname eq $chr && $rnext eq 'chrM')
		{
			if ($pos >= $brk_pt1 && $pos <= $brk_pt2)
				{
				print foutname1 ">$rname\t$pos\t$rnext\t$pnext\n$seq\n";
				}
			}
			
		elsif ($rname eq 'chrM' && $rnext eq $chr)
		{
			if ($pnext >= $brk_pt1 && $pnext <= $brk_pt2)
				{
				print foutname1 ">$rname\t$pos\t$rnext\t$pnext\n$seq\n";
				}
			}
				
    
		}

close (SAM);
	}
close (foutname1);
}



