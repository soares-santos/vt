#!/usr/bin/perl

my $pi = 3.14159266;

while (<>) {

    my ($alpha, $delta) = split(" ", $_);

    $alpha = $alpha / 180 * $pi;
    $delta = $delta / 180 * $pi;

    my $x = cos($alpha) * cos($delta);
    my $y = sin($alpha) * cos($delta);
    my $z = sin($delta);

    print "$x $y $z\n";
}


