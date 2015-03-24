#!/usr/bin/perl -w

# remove-multiples.pl
# Martin Kilbinger
# 2008
# Removes multiple entries (equal first two columns) from input catalogue
# Usage: remove-multiple.pl input > output

if ($#ARGV!=0) {
  print STDERR "Usage: remove-multiples.pl input\n";
  exit 1;
}

die "File $ARGV[0] does not exist" unless -e $ARGV[0];

`sort -bn $ARGV[0] > tmptmp.tmp`;


open(TMP, "tmptmp.tmp") or die "Could not open file tmptmp.tmp: $!";

# First entry
do {
    @last = split(" ", $_ = <TMP>);
} while (/#/ || $#last<=0);          # Ignore comments and one-columner

print;

$nequal = 0;
do {

    $equal = 0;

    do {
	# Get next entry while equal to last
	do {
	    @F = split(" ", $_ = <TMP>);
	} while (/#/ || ($#F<=0 && !eof));

 	# Compare with previous
	if ($last[0]==$F[0] && $last[1]==$F[1]) {
	    $nequal ++;
	    $equal = 1;
	} else {
	    $equal = 0;
	}
    } while ($equal==1 && $_);

    print if $#F>=1;

    # Update last entry
    @last = @F;

    goto count if eof;

} while ($_);


close TMP;

count:
print STDERR "$nequal multiples found and removed (so that Martin's crappy code does not crash)\n"; 

unlink "tmptmp.tmp";

# end
