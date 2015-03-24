#!/usr/bin/perl -w

# stuff.pm
# Martin Kilbinger 2011
# Perl 'include' file
# Include from any perl script as:
#  use lib("$ENV{'HOME'}/perl");
#  use stuff;


# Note: Don't add readdir to Fatal list!
use Fatal qw/ open close opendir closedir mkdir /;
use File::Basename;
use Getopt::Std;
use Cwd;


use warnings FATAL => 'all';


# Constants
$arcmin = 2.9088821e-4;
$arcsec = 4.84813681e-06;



# Runs a command and prints it on stderr
# Exits if command does not return 0
# Default arguments: {verbose => 1, stop => 1, file => undef, run => 1}
sub run_cmd {
    my ($command, $opt_arg) = @_;

    # Default arguments
    $opt_arg->{'verbose'} = 1 unless defined $opt_arg->{'verbose'};
    $opt_arg->{'stop'}    = 1 unless defined $opt_arg->{'stop'};
    $opt_arg->{'run'}     = 1 unless defined $opt_arg->{'run'};

    if (! defined $opt_arg->{'file'} or ! -s $opt_arg->{'file'}) {
        print STDERR "*** Running '$command' ***\n" if $opt_arg->{'verbose'};
        $res = system("$command") if $opt_arg->{'run'};
        $ex  = $?;

        if ($ex!=0) {
            if ($opt_arg->{'verbose'}) {
                my $str = $opt_arg->{'stop'} ? "Stopping" : "Continuing";
                print STDERR "Last command returned error code $ex. $str '" .
                            basename($0) . "'.\n";
            }
            $ex = -1 if $ex == 256;
            exit $ex if $opt_arg->{'stop'};
        }
    } else {
        print STDERR "*** Skipping '$command', non-empty file '$opt_arg->{'file'}' exists ***\n" if $opt_arg->{'verbose'}; 
    }

    #print "run_cmd returning $res\n" if $opt_arg->{'verbose'};
    return $res;
}

sub read_table {
    my ($name) = @_;

    die "read_table: No file name specified" if not defined $name;

    open(my $in_fh, "$name");
    my $i = 0;
    my $ncol = -1;
    my @A = ();
    while (<$in_fh>) {
        next if /#/;
        my @F = split(" ", $_);
        if ($i == 0) {
            $ncol = $#F + 1;
        } else {
            die "Inconsistent file '$name' at line $i ($ncol != ", $#F + 1 , ")" unless $ncol == $#F + 1;
        }
        foreach my $m (0 .. $ncol - 1) {
            $A[$i][$m] = $F[$m];
        }
        $i ++;
    }
    close $in_fh;

    die "Empty file '$name'" if $ncol == -1;

    my $nlines = $i;

    return \@A, $nlines, $ncol;
}

# Returns a hash from a list with indices corresponding to tags.
# Use as:
#  my %rec = %{get_tags_idx(@header)}
sub get_tags_idx {

  my @header = @_;

  my %rec = ();
  foreach my $i (0 ..$#header) {
    my $tag = defined $rec{$header[$i]} ? $header[$i] . "_$i" : $header[$i];
    $rec{$tag} = $i;
  }

  return \%rec;
}

# Returns all lines from a file and the number of non-comment lines
sub read_file {
  my ($fname) = @_;

  open(my $in_fh, "$fname");
  my @lines = ();
  my $i = 0;
  while (<$in_fh>) {
    next if /^#/;
    $lines[$i] = $_;
    $i ++;
  }

  return @lines;
}


# Returns all file names in directory $dir which match $name_pattern
sub read_files {
    my ($name_pattern, $dir, $opt_arg) = @_;
    $dir = "." unless defined $dir;

    $opt_arg->{'verbose'} = defined $opt_arg->{'verbose'} ? $opt_arg->{'verbose'} : 0;

    print STDERR "Looking for files '$name_pattern' in dir '$dir'\n" if $opt_arg->{'verbose'};

    opendir(my $dir_fh, "$dir");
    my @FILES = grep {/^$name_pattern$/} readdir($dir_fh);
    closedir $dir_fh;

    return @FILES;
}

# Returns the names of all subdirectories.
# New: argument list as parent directories
# (default ".")
sub read_subdirs {
    my @DIR = @_;

    my @subdirs = ();
    # Parent directories
    @DIR = (".") if $#DIR == -1;

    my $dir_fh = undef;
    foreach $dir (@DIR) {
        opendir($dir_fh, $dir);
        # Alle Unterverzeichnisse
        my @to_add = grep {-d "$dir/$_" && ! /^\.{1,2}$/} readdir($dir_fh);
        foreach $t (0 .. $#to_add) {
            $to_add[$t] = "$dir/$to_add[$t]"; 
        }
        #print "$dir: @to_add\n";
        closedir($dir_fh);
        @subdirs = (@subdirs, @to_add);
    }

    @subdirs = sort(@subdirs);
    return @subdirs;
}

# Write an 'ahtena' config file ('config_tree')
sub write_config_athena {
    my ($opt_arg) = @_;

    #print "$opt_arg\n";
    #my(%globals) = %{$opt_arg}; while ( ($key, $value) = each %globals) { print "$key => $value\n"; }

    my $out_name  = defined $opt_arg->{'out_name'} ? $opt_arg->{'out_name'} : "config_tree";
    my $cat1      = defined $opt_arg->{'GALCAT1'}  ? $opt_arg->{'GALCAT1'}  : "gal";
    my $cat2      = defined $opt_arg->{'GALCAT2'}  ? $opt_arg->{'GALCAT2'}  : "-";
    my $wcorr     = defined $opt_arg->{'WCORR'}    ? $opt_arg->{'WCORR'}    : "1";
    my $format    = defined $opt_arg->{'SFORMAT'}  ? $opt_arg->{'SFORMAT'}  : "standard";
    my $sc_in     = defined $opt_arg->{'SCOORD_INPUT'} ? $opt_arg->{'SCOORD_INPUT'} : "arcmin";
    my $sc_out    = defined $opt_arg->{'SCOORD_OUTPUT'} ? $opt_arg->{'SCOORD_OUTPUT'} : "arcmin";
    my $thmin     = defined $opt_arg->{'THMIN'}    ? $opt_arg->{'THMIN'}    : 0.05;
    my $thmax     = defined $opt_arg->{'THMAX'}    ? $opt_arg->{'THMAX'}    : 240.0;
    my $nth       = defined $opt_arg->{'NTH'}      ? $opt_arg->{'NTH'}      : 1000;
    my $bintype   = defined $opt_arg->{'BINTYPE'}  ? $opt_arg->{'BINTYPE'}  : "LIN";
    my $radec     = defined $opt_arg->{'RADEC'}    ? $opt_arg->{'RADEC'}    : 0;
    my $oath      = defined $opt_arg->{'OATH'}     ? $opt_arg->{'OATH'}     : 0.03;
    my $serror    = defined $opt_arg->{'SERROR'}   ? $opt_arg->{'SERROR'}   : "none";
    my $nresample = defined $opt_arg->{'NRESAMPLE'} ? $opt_arg->{'NRESAMPLE'} : "0 0";

    open(my $fh, ">$out_name");
    print "Writing $out_name\n";
    print {$fh} "### -*- sh -*- ###\n";
    print {$fh} "### Config file for tree code 'athena'\n\n";
    print {$fh} "GALCAT1\t\t$cat1\n";
    print {$fh} "GALCAT2\t\t$cat2\n";
    print {$fh} "WCORR\t\t$wcorr\n";
    print {$fh} "SFORMAT\t\t$format\n";
    print {$fh} "SCOORD_INPUT\t$sc_in\n";
    print {$fh} "SCOORD_OUTPUT\t$sc_out\n";
    print {$fh} "THMIN\t\t$thmin\n";
    print {$fh} "THMAX\t\t$thmax\n";
    print {$fh} "NTH\t\t$nth\n";
    print {$fh} "BINTYPE\t\t$bintype\n";
    print {$fh} "RADEC\t\t$radec\n";
    print {$fh} "OATH\t\t$oath\n";
    print {$fh} "SERROR\t\t$serror\n";
    print {$fh} "NRESAMPLE\t$nresample\n" unless $serror eq "none";
    close $fh;
}

sub write_log {
    my $name = $_[0];
    shift @_;
    my @ARGV = @_;

    open(my $fh, ">$name");
    print {$fh} "$0";
    foreach $argv (@ARGV) {
        @largv = split(" ", $argv);
        if ($#largv == 0) {
            print {$fh} " $argv";
        } else {
            print {$fh} " \"$argv\"";
        }
    }
    print {$fh} "\n";
    close $fh;
}

sub get_cmd_line {
    my @ARGV = @_;

    my $cmd = "$0";
    foreach $argv (@ARGV) {

        # Check whether argument is a white-space-separated list
        @largv = split(" ", $argv);
        if ($#largv == 0) {
            $cmd = $cmd . " $argv";
        } else {
            $cmd = $cmd . " \"$argv\"";
        }
    }
    $cmd = $cmd . "\n";

    return $cmd; 
}

sub write_cmd_line {
    my ($name, $cmd) = @_;
    open(my $fh, ">$name");
    print {$fh} "$cmd";
    close $fh;
}

# Also defined in POSIX
sub log10 {
    my ($x) = @_;
    return log($x) / log(10);
}

# From the Perl Cookbook (http://www.oreilly.com/catalog/perlckbk2)
sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}

# Color support
use Term::ANSIColor;
#print "nok" if not defined color;
#color('red');


# Return value for module
1;
