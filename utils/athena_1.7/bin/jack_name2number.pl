#!/usr/bin/perl -w 

# jack_name2number.pl
# Martin Kilbinger 2012
# Transforms catalogue from Jackknife sample name
# to Jackknife sample number format


use stuff;

my %options = ();
getopts("o:c:h", \%options);

usage(0) if defined $options{h};
usage(-1) if $#ARGV != 0;

my $cat_name = $ARGV[0];
my $col      = defined $options{c} ? $options{c} : 3;

# Go through file, first round
my %name = ();
my $jk_num = 0;
open(my $in_fh, "$cat_name");
while (<$in_fh>) {
    next if /#/;
    my @F = split(" ", $_);


    # Set Jackknife number for each name
    if (! defined $name{$F[$col]}) {
        $name{$F[$col]} = $jk_num;
        $jk_num ++;
    }
}
close $in_fh;

#while ( ($key, $value) = each %name ) { print "$key => $value\n"; } exit 0;


# Write output, replace name by number
my $out_fh = *STDOUT;

open($out_fh, ">$options{o}") if defined $options{o};
open($in_fh, "$cat_name");
while (<$in_fh>) {
    (print && next) if /#/;
    my @F = split(" ", $_);

    foreach my $i (0 .. $#F) {
        if ($i == $col) { print {$out_fh} "$name{$F[$col]} "; }
        else { print {$out_fh} "$F[$i] "; }
    }
    print {$out_fh} "\n";
}
close $in_fh;
close $out_fh;


sub usage {
    my ($ex) = @_;

    print STDERR "Usage: jack_name2number.pl [OPTIONS] CAT\n";
    print STDERR "OPTIONS\n";
    print STDERR "  CAT         Input catalogue in JK sample name format\n";
    print STDERR "  -c COL      Name in column COL (default 3, corresponding to 'position_jack_num' format)\n";
    print STDERR "  -o OUT      Output file OUT (default stdout)\n";
    print STDERR "  -h          This message\n";


    exit $ex if defined $ex;
}

