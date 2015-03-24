#!/usr/bin/perl -w

if ($#ARGV!=1) {
    print STDERR "Usage: jean2cat.pl in out\n";
    exit 1;
}

$c_x    = 3;
$c_y    = 4;
$c_flag = 6;
$c_i    = 10;

open(IN, "$ARGV[0]") or die "$!";
open(OUT, ">$ARGV[1]") or die "$!";

while (<IN>) {
    next if /#/;

    @F = split(" ", $_);

    # Bad/masked object
    next if $F[$c_flag]!=0;

    # Magnitude cuts
    next if $F[$c_i] >= 22.5;
    #next if $F[$c_i] > 23;

    print OUT "$F[$c_x] $F[$c_y]\n";

}
close IN;
close OUT;

# end


