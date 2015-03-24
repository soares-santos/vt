#!/usr/bin/perl -w 

# jack2jack.pl
# Martin Kilbinger 2012
# Transforms catalogues between different jackknife formats
# (Jackknife sample list; Jackknife sample number)

use FindBin;
use lib $FindBin::Bin;
use stuff;

my %options = ();
getopts("h", \%options);

usage(0) if defined $options{h};
usage(-1) if $#ARGV != 0;

my $cat_name = $ARGV[0];


# Go through file, first round
my $ncol = 0;
my $njmax = -1;
open(my $in_fh, "$cat_name");
while (<$in_fh>) {
    next if /#/;
    my @F = split(" ", $_);

    $ncol = $#F if $ncol == 0;
    last if $ncol != 3;  # Catalogue is in Jackknife sample list format -> end scanning

    # Here: Catalogue is in Jackknife sample number format -> scan for largest number
    my $nj_num = $F[3];
    $njmax = $nj_num if $njmax < $nj_num;
}
close $in_fh;



my $njack = $ncol - 2;   # For JK sample list mode only

print STDERR "njmax = $njmax\n";
print STDERR "njack = $njack\n";
print STDERR "ncol  = $ncol\n";

# Reopen file
open($in_fh, "$cat_name");
while (<$in_fh>) {
    next if /#/;
    my @F = split(" ", $_);

    print "$F[0] $F[1] $F[2]";

    if ($ncol == 3) {

	    # Input catalogue is Jackknife sample number format
	    my $found = 0;
	    foreach my $i (0 .. $njmax) {
	        if ($i == $F[3]) {
		        print " 0";
		        $found = 1;
	        } else {
		        print " 1";
	        }
	    }
	    if ($found == 0) {
	        print STDERR "not found\n";
	        out_line(@F);
	        exit 1;
	    }

    } else {

	    # Input catalogue is Jackknife sample list format
	    my $found = 0;
	    foreach my $i (0 .. $njack - 1) {
	        if ($F[$i + 3] == 0) {
		        print " $i";
		        $found ++;
	        }
	    }
	    if ($found != 1) {
	        print STDERR "Not found ($found)\n";
	        out_line(@F);
	        exit 2;
	    }

    }

    print "\n";

}
close $in_fh;


sub out_line {
    my @F = @_;

    foreach my $f (@F) {
	print STDERR "$f ";
    }
    print STDERR "\n";
}

sub usage {
    my ($ex) = @_;

    print STDERR "Usage: jack2jack.pl [OPTIONS] CAT\n";
    print STDERR "OPTIONS\n";
    print STDERR "  CAT         Input catalogue in JK sample list or number format\n";
    print STDERR "               (The correct format is determined by the script.)\n";
    print STDERR "  -h          This message\n";


    exit $ex if defined $ex;
}

