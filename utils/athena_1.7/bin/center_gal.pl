#!/usr/bin/perl -w

# center_gal.pl v1.1
# Martin Kilbinger 2007
# Calculates the center of all galaxies from a catalog using ra and dec.

use Math::Trig;
use Fatal qw/ open /;
use Getopt::Long;
%options=();

GetOptions("ra=i"     =>\$options{ra},
   	       "dec=i"    =>\$options{dec},
           "head!"    =>\$options{head});

$c_ra  = defined $options{ra} ? $options{ra}  : 0;
$c_dec = defined $options{dec}? $options{dec} : 1;
$options{head} = 1 if !defined $options{head};


if ($#ARGV!=0) {
    printf STDERR "Usage: center_gal.pl [-ra alpha_col -dec delta_col] [-[no]head] catalog\n";
    exit(5);
}

die "File '$ARGV[0]' not found" unless -e $ARGV[0];
die "File '$ARGV[0]' empty" unless -s $ARGV[0];


if ($options{head}) {
    print "# alpha delta x y z\n";
    print "# center of all galaxies\n";
}

@r = (0,0,0);

$n = 0;
while (<>) {

    @F = split(" ", $_);
    next if /#/;
    next if $#F==0;

    ($alpha, $delta) = ($F[$c_ra], $F[$c_dec]);

    $cd    = cosd($delta);
    $r[0] += cosd($alpha)*$cd;
    $r[1] += sind($alpha)*$cd;
    $r[2] += sind($delta);
    $n++;
}

for $k (0, 1, 2) {
    $r[$k] = $r[$k]/$n;
}

$absr = sqrt($r[0]**2 + $r[1]**2 + $r[2]**2);
for $k (0, 1, 2) {
    $r[$k] = $r[$k]/$absr;
}

$alpham = atan2($r[1], $r[0])*180/pi;
if ($alpham<0) { $alpham += 360; }
$deltam = asin($r[2])*180/pi;

printf "%.4f %.4f %f %f %f\n", $alpham, $deltam, $r[0], $r[1], $r[2];


# trig functions for arguments in deg
sub cosd { cos($_[0]/180.0*pi); }
sub sind { sin($_[0]/180.0*pi); }
