#!/usr/bin/perl -w

# cat2gal.pl
#
# Martin Kilbinger 2007
# v1.1 (2008)
# v1.2 (2012): Jackknife number added
#
# Reads a catalog and writes a file in a format suitable
# for the tree code 'athena'.
# Writes
# x y [e1 e2] w [jk].
# Also performs coordinate projections.
# For the rotations see Calabretta&Greisen (2002).
#
# This script is part of the athena package,
# www.cosmostat.org/athena

use Math::Trig;
use Fatal qw/ open /;
use List::Util qw/ max /;

use Getopt::Long;
%options=();

GetOptions("x=i"        => \$options{x},
	   "y=i"            => \$options{y},
	   "coord_input=s"  => \$options{coord_input},
	   "coord_output=s" => \$options{coord_output},
	   "e1=i"           => \$options{e1},
	   "e2=i"           => \$options{e2},
	   "w=i"            => \$options{w},
       "jk=i"           => \$options{jk},
	   "alphac=f"       => \$options{alphac},
	   "deltac=f"       => \$options{deltac},
	   "q"              => \$options{q},
	   "h"              => \$options{h}
);


$coord_input  = defined $options{coord_input} ? $options{coord_input} : "rad";
$coord_output = defined $options{coord_output} ? $options{coord_output} : "rad";

$c_x       = defined $options{x}  ? $options{x}  :  0;
$c_y       = defined $options{y}  ? $options{y}  :  1;
$c_e1      = defined $options{e1} ? $options{e1} : -1;
$c_e2      = defined $options{e2} ? $options{e2} : -1;
$c_w       = defined $options{w}  ? $options{w}  : -1;
$c_jk      = defined $options{jk} ? $options{jk} : -1;
$c_max     = max($c_x, $c_y, $c_e1, $c_e2, $c_w, $c_jk);
$quiet     = defined $options{q}  ? 1            :  0;


usage() if defined $options{h};
usage() if ($#ARGV!=1);

die "Only one of -e1, -e2 given" if ($c_e1<0 && $c_e2>=0) || ($c_e2<0 && $c_e1>=0);

$cat  = $ARGV[0];
$mode = $ARGV[1];


# Conversion factors
$conv{rad}    = 1;
$conv{arcsec} = 4.84813681e-06;
$conv{arcmin} = 2.90888209e-04;
$conv{deg}    = 0.01745329;

die "Coordinate '$coord_input' not defined" unless defined $conv{$coord_input};
die "Coordinate '$coord_output' not defined" unless defined $conv{$coord_output};


if ($mode eq "tan") {
  print STDERR "tan (gnonomic) projection\n" unless $quiet;
} elsif ($mode eq "no") {
  print STDERR "no projection\n" unless $quiet;
} elsif ($mode eq "skynew" || $mode eq "deltac") {
  print STDERR "cos(delta_fieldcenter) projection\n" unless $quiet;
} else {
  print STDERR "wrong mode\n";
  exit(2);
}

# Determine number of galaxies
open (CAT, "$cat") or die "could not open file $cat: $!";

# Default values if not in catalogue
$erot1 = $erot2 = 0.0;
$w = 1.0;

# Calculate field center if necessary
$calc_center = 0;
if ($mode eq "skynew" || $mode eq "deltac") {

  if (! defined $options{deltac}) {
    $calc_center = 1;
  } else {
    $deltac = $options{deltac};
  }

} elsif ($mode eq "tan") {

  if (! defined $options{alphac} && ! defined $options{deltac}) {
    $calc_center = 1;
  } else {
    $alphac = $options{alphac}*$conv{$coord_input};
    $deltac = $options{deltac}*$conv{$coord_input};
    print STDERR "Fieldcenter (from command line) = ($alphac, $deltac)\n" unless $quiet;
  }

}

if ($calc_center) {
  $F = `center_gal.pl -ra $c_x -dec $c_y -nohead $cat`;
  ($alphac, $deltac, $xc, $yc, $zc) = split(" ", $F);
  $alphac = $alphac*$conv{$coord_input};
  $deltac = $deltac*$conv{$coord_input};
  print STDERR "Fieldcenter (calculated) = ($alphac, $deltac)\n" unless $quiet;
}


# Transform coordinates
open(CAT, "$cat");
my ($xmin, $xmax, $ymin, $ymax) = (1e30, -1e30, 1e30, -1e30);
while (<CAT>) {

	@F = split(" ", $_);
    next if /#/;			     	# Header
    (print && next) if $#F==0;		# Number of galaxies (optional)

    die "Column #$c_max not available in input catalogue '$cat'" if $#F < $c_max;

	$x       = $F[$c_x];
	$y       = $F[$c_y];
	$erot1   = $F[$c_e1] if $c_e1>=0;
	$erot2   = $F[$c_e2] if $c_e2>=0;
	$w       = $F[$c_w]  if defined $options{w};
    $jk      = $F[$c_jk] if defined $options{jk};


	# Input coord -> rad
	$x = $x*$conv{$coord_input};
	$y = $y*$conv{$coord_input};


	if ($mode eq "tan") {

	  $alpha = $x;
	  $delta = $y;
	  ($x, $y, $z)     = ad2xyz($alpha, $delta);

	  # Rotate center of field (alphac, deltac) to (0, 0)
	  ($x, $y, $z)     = rotate(-$alphac, 2, $x, $y, $z);
	  ($x, $y, $z)     = rotate(-$deltac, 1, $x, $y, $z);

	  ($alpha, $delta) = xyz2ad($x, $y, $z);

	  # Gnomonic projection: r'/r = tan(phi)/sin(phi)
	  $cosphi = cos($alpha)*cos($delta);
	  $x      = $alpha/$cosphi;
	  $y      = $delta/$cosphi;

	} elsif ($mode eq "skynew" || $mode eq "deltac") {

	  # Multiply alpha with the cosine of the field center
	  $x = $x*cos($deltac*$conv{$coord_input});

	}

    # Store min and max positions
    $xmin = $x if $xmin > $x;
    $ymin = $y if $ymin > $y;
    $xmax = $x if $xmax < $x;
    $ymax = $x if $ymax < $x;

	# Rad -> output coord
	$x = $x/$conv{$coord_output};
	$y = $y/$conv{$coord_output};


    printf "% .8f % .8f", $x, $y;
    printf " % .8f % .8f", $erot1, $erot2 if $c_e1>0 and $c_e2>0;
    printf " % .8f", $w;
    printf " %4d", $jk if $c_jk>0;
    printf "\n";

}
close CAT;

# Check whether maximum Cartesian distance is not > 180 deg
if (! ($mode eq "no") ) {
    my $pi = 3.14159266;
    my $dmaxsqr = ($xmax - $xmin)**2 + ($ymax - $ymin)**2;
    if ($dmaxsqr > $pi*$pi) {
        #print "$xmax $xmin $ymax $ymin\n";
        my $dmax = sqrt($dmaxsqr);
        my $dmax_deg = $dmax / $pi * 180;
        die "Maximum distance (", $dmax, " = ", $dmax_deg, " deg) is larger than 180 deg, projection to Cartesian coordinates is not possible"
            if $dmaxsqr > $pi*$pi;
    }
}

# Avoid 'unused' warnings
$alphac = $xc = $yc = $zc = 0.0;


sub xyz2ad {

  my ($x, $y, $z) = @_;

  my $delta = asin($z);
  my $alpha = atan2($y, $x);

  return ($alpha, $delta);
}

sub ad2xyz {
  my ($alpha, $delta) = @_;

  my $x = cos($alpha)*cos($delta);
  my $y = sin($alpha)*cos($delta);
  my $z = sin($delta);

  return ($x, $y, $z);
}

sub rotate {
  my ($phi, $axis, $x, $y, $z) = @_;

  my $cp = cos($phi);
  my $sp = sin($phi);

  if ($axis==2) {
    $xn = $x*$cp - $y*$sp;
    $yn = $x*$sp + $y*$cp;
    $zn = $z;
  } elsif ($axis==1) {
    $xn = $x*$cp - $z*$sp;
    $yn = $y;
    $zn = $x*$sp + $z*$cp;
  } elsif ($axis==0) {
    $xn = $x;
    $yn = $y*$cp - $z*$sp;
    $zn = $y*$sp + $z*$cp;
  }

  return ($xn, $yn, $zn);
}


sub usage {
    print STDERR "Usage: cat2gal.pl [OPTIONS] CAT TYPE\n";
    print STDERR "OPTIONS:\n";
    print STDERR "  -coord_input COORD        Input coordinates COORD = {'rad' (default), 'arcsec', 'arcmin', 'deg'}\n";
    print STDERR "  -coord_output COORD       Output coordinates COORD = {'rad' (default), 'arcsec', 'arcmin', 'deg'}\n";
    print STDERR "  -x COL                    Column COL of x-coordinate (default: 0)\n";
    print STDERR "  -y COL                    Column COL of y-coordinate (default: 1)\n";
    print STDERR "  -e1 COL                   Column COL of ellipticity, first component (default: unused)\n";
    print STDERR "  -e2 COL                   Column COL of ellipticity, second component (default: unused)\n";
    print STDERR "  -w COL                    Column COL of weight (default: unused, all weights set to 1)\n";
    print STDERR "  -jk COL                   Column COL of jackknife number (default: unused)\n";
    print STDERR "  -alphac ALPHAC            Right ascension of field center ALPHAC (for 'tan' projection)\n";
    print STDERR "  -deltac DELTAC            Declination of field center DELTAC (for 'deltac' and 'tan' projection)\n";
    print STDERR "  -q                        Quiet(er) mode\n";
    print STDERR "  CAT                       Input catalogue CAT\n";
    print STDERR "  TYPE                      Projection type TYPE:\n";
    print STDERR "                              'no':      No projection, keep input (x,y) unchanged\n";
    print STDERR "                              'tan':     Gnomonic projection\n";
    print STDERR "                              'deltac':  'cos(delta_c)'-projection\n";
    print STDERR "\nIf TPYE='tan' or 'deltac':\n";
    print STDERR "   Input coordinates (x,y) are assumed to be spherical coordinates (ra, dec).\n";
    print STDERR "   If DELTAC (and ALPHAC) are not given, they are calculated by calling the script 'center_gal.pl.\n";
    print STDERR "   If weight column (option '-w') not given, a weight of unity is written\n";
    exit(4);
}

# endf
