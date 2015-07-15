#!/bin/perl

my $INFILE = $ARGV[0];
my $WINDOW = $ARGV[1] || 100;
my $MINCOUNTS = $ARGV[2] || 0;
my $LABEL = $ARGV[3] || "Label";

# perl ../../../scripts/C_context_window.pl CHG_context_Bismark_se.1321246134.txt 10000 0 CHG_B73

# HWI-ST261:400:B0A1EABXX:8:2208:11321:57205 1:N:0:ATCACG	-	7	16761635	1	x
# HWI-ST261:400:B0A1EABXX:8:2103:16659:181377 1:N:0:ATCACG	+	1	248757208	X

my %dat;
my $OUTFILE = "$INFILE.$WINDOW.bed";
# Read in data
open (IN, $INFILE) or die;
open (OUT, ">$OUTFILE") or die;

print STDERR "Reading in report file...\n";
my $ts = localtime();
print STDERR $ts, "\n";

while (my $l = <IN>) {

	chomp($l);
	my @f = split(/\t/, $l);
	
	my $chr = $f[0];
	my $pos = $f[1];
	
	my $bin = int($pos / $WINDOW);
	
	$dat{$chr}->{$bin}->{'total'}++;
	
	#print STDERR "$chr\t$bin\t$met\n";
	
}

print STDERR "Printing out window file...\n";
$ts = localtime();
print STDERR $ts, "\n";


# Emit window bigBed
print STDERR "Printout time\n";
for my $C (sort {$a <=> $b} keys %dat) {
	
	for my $W (sort {$a <=> $b} keys %{ $dat{$C} }) {
	#while(my ($W, $value) = each(%{ $dat{$C} })) {
			
		# Don't compute over windows with too few C's
		if ($dat{$C}->{$W}->{'total'} > $MINCOUNTS) {

	print OUT $C, "\t", ($W * $WINDOW), "\t", (($W + 1) * $WINDOW) - 1, "\t", $dat{$C}->{$W}->{'total'}, "\n";
			}
		}
	}

close OUT;
close IN;
