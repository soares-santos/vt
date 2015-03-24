#!/usr/bin/perl -w

# plot_cat_jack_list.pl (part of athena)
# Martin Kilbinger
# 2012
# Plots galaxies in different colors, encoding their Jackknife
# sample number. Input is a catalogue in Jackknife sample list
# format.


use FindBin;
use lib $FindBin::Bin;
use stuff;

my %options = ();
getopts("e:kh", \%options);

usage(0) if defined $options{h};
usage(-1) if $#ARGV != 0;

my $ev       = defined $options{e} ? $options{e} : 2;
my $cat_name = $ARGV[0];
my $out_base = "$cat_name";

my $ncol  = 0;
open(my $tmp_fh, "$cat_name");
while (<$tmp_fh>) {
    next if /#/;
    my @F = split(" ", $_);
    $ncol = $#F;
    last;
}
close $tmp_fh;

my $njack = $ncol - 2;

if ($njack == 1) {
    print STDERR "Warning: Only one Jackknife column found. Catalogue is maybe in Jackknife sample number format.\n";
    print STDERR "Use 'jack2jack.pl $ARGV[0] > $ARGV[0].jack_list' to create a catalogue in Jackknife list format.\n";
} else {
    print STDERR "$njack Jackknife samples found\n";
}

open(my $out_fh, ">$out_base.gnu");

print {$out_fh} "set terminal post eps enhanced color 'Times-Roman' 20\n";
print {$out_fh} "set output '$out_base.eps'\n";
print {$out_fh} "unset key \n";
print {$out_fh} "pl '$cat_name' ev $ev ps 1 lc 0";
foreach my $j (4 .. $ncol + 1) {
    print {$out_fh} ", '' ev $ev u 1:(\$$j == 0 ? \$2 : 1/0) ps 1";
}
print {$out_fh} "\n";
print {$out_fh} "set output\n";
close $out_fh;

run_cmd("gnuplot $out_base.gnu");

unlink "$out_base.gnu" unless defined $options{k};


sub usage {
    my ($ex) = @_;

    print STDERR "Usage: plot_cat_jack.pl [OPTIONS] CAT\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  CAT         Catalogue in Jackknife sample list format (JK indices in columns 2 .. NJ+2)\n";
    print STDERR "  -e EV       Plot every EV point\n";
    print STDERR "  -k          Don't delete gnuplot script\n";
    print STDERR "  -h          This message\n";

    exit $ex if defined $ex;
}
