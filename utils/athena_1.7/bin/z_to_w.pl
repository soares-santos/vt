#!/usr/bin/perl -w

my $cat_in_name = $ARGV[0];

open(my $out_fh, ">comovtmp.i");
print {$out_fh} "include, \"stuff.i\"\n";
print {$out_fh} "include, \"distance.i\"\n";

print {$out_fh} "Omega_m = 0.7\n";
print {$out_fh} "Omega_L = 1 - Omega_m\n";
print {$out_fh} "h = 0.7\n";
print {$out_fh} "H_0 = 100 * h\n";

print {$out_fh} "read_table, \"$cat_in_name\", cat\n";
print {$out_fh} "x = cat(1,)\n";
print {$out_fh} "y = cat(2,)\n";
print {$out_fh} "z = cat(3,)\n";

# deg -> rad
print {$out_fh} "x = x / 180 * pi\n";
print {$out_fh} "y = y / 180 * pi\n";

# redshift -> com. distance
print {$out_fh} "Z = w(z)\n";
print {$out_fh} "Zm = Z(avg)\n";

# rad -> com. distance (X = alpha * Z)
print {$out_fh} "X = (x - x(avg)) * Z\n";
print {$out_fh} "Y = (y - y(avg)) * Z\n";

# Leave angular coordinates
#print {$out_fh} "X = x\n";
#print {$out_fh} "Y = y\n";

my $cat_out_name = "$cat_in_name" . ".comoving";

print {$out_fh} "f = open(\"$cat_out_name\", \"w\")\n";
print {$out_fh} "write, f, format=\"% f % f % f 1.0\\n\", X, Y, Z\n";
print {$out_fh} "close, f\n";
print {$out_fh} "quit\n";

close $out_fh;

system("yorick -batch comovtmp.i");

# end
