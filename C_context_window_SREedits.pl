#!/bin/perl

#Windowing script for Bismark methylation extraction data
#provide bismark file, window size, coverage minimum (0 = no minimum) and track title

my $INFILE = $ARGV[0];
my $WINDOW = $ARGV[1] || 100;
my $MINCOUNTS = $ARGV[2] || 0;
my $LABEL = $ARGV[3] || "Label";

# perl ../../../scripts/C_context_window.pl CHG_context_Bismark_se.1321246134.txt 10000 0 CHG_B73


#bismark methylation extractor output calls methylation by the '+' or '-' flag in each file as described below:
# HWI-ST261:400:B0A1EABXX:8:2208:11321:57205 1:N:0:ATCACG	-	7	16761635	1	x
# HWI-ST261:400:B0A1EABXX:8:2103:16659:181377 1:N:0:ATCACG	+	1	248757208	X

my %dat;
my $OUTFILE = "$INFILE.$WINDOW.wig";
# Read in data
open (IN, $INFILE) or die;
open (OUT, ">$OUTFILE") or die;

print STDERR "Reading in report file...\n";
my $ts = localtime();
print STDERR $ts, "\n";

#go through each line in the file
while (my $l = <IN>) {

	chomp($l);
	my @f = split(/\t/, $l);
	#grab the chromosome, position, and methylation call
	my $chr = $f[2];
	my $pos = $f[3];
	my $met = $f[1];
	
	#define the 'bin' by dividing position by windowsize and making it an integer
	my $bin = int($pos / $WINDOW);
	my $sites = $pos;
	
	#store your data as nested hashes
	$dat{$chr}->{$bin}->{'total'}++;
	#the positions are stored as well to allow identification of the number of sites that are being looked at
	$dat{$chr}->{$bin}->{$sites}->{'sitecount'}++;

	if ($met eq '+') {
		$dat{$chr}->{$bin}->{'met'}++;		
	}	
}

	

print STDERR "Printing out window file...\n";
$ts = localtime();
print STDERR $ts, "\n";

print OUT "track type=bedGraph name=", $LABEL, "_", $WINDOW, "_bp\n";

# Emit window bigBed
print STDERR "Printout time\n";
#sort your hash of hashes to take the first chromosome and first bin...
for my $C (sort {$a <=> $b} keys %dat) {
	for my $W (sort {$a <=> $b} keys %{ $dat{$C} })  {
			
		# Don't compute over windows with too few C's
		if ($dat{$C}->{$W}->{'total'} > $MINCOUNTS) {
		
			#define the fraction of methylation by taking the number of methylation reads / total reads in your window
			my $frac = int( 100 * ($dat{$C}->{$W}->{'met'} / $dat{$C}->{$W}->{'total'}));
			
			#if there is methylated reads in the window
			if ($C ne '' && $dat{$C}->{$W}->{'met'}) {
					#print out everything. To get the number of sites looked at, I am taking the number of keys in each bin minus two. This is to remove the 'total' and 'met' keys
							print OUT $C, "\t", ($W * $WINDOW), "\t", (($W + 1) * $WINDOW) - 1, "\t", $frac, "\t", $dat{$C}->{$W}->{'met'}, "\t", $dat{$C}->{$W}->{'total'}-$dat{$C}->{$W}->{'met'}, "\t", $dat{$C}->{$W}->{'total'}, "\t", scalar(keys $dat{$C}->{$W})-2, "\n";
			} else {if ($C ne ''){
				#same as above, however how we only subtract one 'total' as there is no key for 'met' if there were not reads for it
				print OUT $C, "\t",($W * $WINDOW), "\t", (($W + 1) * $WINDOW) - 1, "\t", $frac, "\t", 0, "\t", $dat{$C}->{$W}->{'total'}-$dat{$C}->{$W}->{'met'}, "\t", $dat{$C}->{$W}->{'total'}, "\t", scalar(keys $dat{$C}->{$W})-1, "\n";
			}}
		}
	}

}

close OUT;
close IN;
