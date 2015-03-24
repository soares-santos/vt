#!/usr/bin/perl -w


##############################################################
# Martin Kilbinger 2008-2013                                 #
# Calculates the spatial correlation function of galaxies.   #
# Calls J. Coupon's venice to create a random catalogue and  #
# then the tree code to calculate the correlation functions. #
##############################################################


# The version has to match the one in tree/main.h
$version = 1.7;


use FindBin;
use Class::Struct;
use List::Util qw/ min max /;
use Fatal qw/ open rename /;
use Getopt::Std;
use Thread qw/ async /;
use File::Basename;


# Path to athena/bin. This script is $athena_bin/woftheta_xcorr.pl
$athena_bin = $FindBin::Bin;

# Check version consistency
my $v = qx($athena_bin/athena -v);
chomp($v);
die "athena version ($v) inconsistent with woftheta_xcorr.pl version ($version)" unless $version eq $v;


##############################
### Command line arguments ###
##############################

print STDERR "woftheta_xcorr.pl v$version (author: Martin Kilbinger)\n";
%options = ();
getopts("c:anp:dRrqh", \%options);

usage(0) if defined $options{h};
usage(1) if $#ARGV!=-1;

$nproc = defined $options{p} ? $options{p} : 1;
$quiet = defined $options{q} ? 1 : 0;
$qflag = $quiet==1 ? "-q" : "";


###################
### Config file ###
###################

# Necessary config file entries (all are scalars except "specifiers" which is a vector)
@conf_entry = (
    "catbase", "catbase2", "path", "path_ran", "rand", "rand2", "mask", "nrand", "nrand2",
    "coord_type", "coord_units", "project", "c_x", "c_y", "specifiers",
    "coord_output", "thmin", "thmax", "nth", "bintype", "oath", "error", "nresample",
);

# Optional config file entries
@conf_entry_opt = (
   "c_jk",
);

my $config_name = defined $options{c} ? $options{c} : "config_woftheta";
read_config($config_name);
check_config();
#out_config();


#######################
### Catalogue names ###
#######################

# Data catalogues
$data1 = "$conf{path}/$conf{catbase}";

# Random catalogue (to be created)
$random1 = "$conf{path_ran}/$conf{rand}";


if (!($conf{catbase2} eq "-")) {

    # Correlations between two catalogue sets (auto- for equal bins, cross- for different bins)
    $ncat = 2;
    $data2   = "$conf{path}/$conf{catbase2}";
    die "Second random catalogue name missing" if $conf{rand2} eq "-";
    die "Second random catalogue name ('$conf{rand2}') has to be distinct from first" if $conf{rand2} eq $conf{rand};
    #die "Number for second random catalogue ('config:nrand2') not valid" if $conf{nrand2} eq "-";
    $random2 = "$conf{path_ran}/$conf{rand2}";

    die "Second random catalogue name ('config:rand2') has to be specified if there are two data catalogues"
      if $conf{rand2} eq "-";
    die "First and second data catalogue cannot have identical names"
      if $conf{catbase} eq $conf{catbase2};
    die "First and second random catalogue cannot have identical names. If you only have one random catalogue, dupliate the file and rename it."
      if $conf{rand} eq $conf{rand2};

} else {

    # Cross-correlations between two catalogue sets
    $ncat    = 1;
    $data2   = "$data1";
    $random2 = "$random1";

    die "Second random catalogue ('config:rand2') given but no second data catalogue ('config:catbase2')"
      unless $conf{rand2} eq "-";

}

# File extensions
$end_final  = $conf{project} eq "none" ? "-orig.dat" : "-proj.dat";

# Catalogue specifiers (redshift ranges until v1.51)
@str = $#specifiers == -1 ? ("") : @specifiers;

### Clean up? ###
if (! defined $options{n}) {
  clean_up($options{a});
}



############################
### Deal with catalogues ###
############################

# Format and transform input data cat to 2pcf input format
$alpha_c_0_data1 = undef;
$delta_c_0_data1 = undef;
$alpha_c_0_data2 = undef if $ncat == 2;
$delta_c_0_data2 = undef if $ncat == 2;
foreach $i (0 .. $#str) {

    # Projection. Calculate center the first time, then use this value for successive bins
    ($atmp, $dtmp) = create_final_cat("$data1$str[$i]", "$end_final", $conf{c_x}, $conf{c_y}, $alpha_c_0_data1, $delta_c_0_data1, $conf{c_jk});
    $alpha_c_0_data1 = $atmp;
    $delta_c_0_data1 = $dtmp;
    # MKDEBUG: The last specifier centre is used for the random catalogue projection. Is this always a good idea?

    if ($ncat == 2) {
        ($atmp, $dtmp) = create_final_cat("$data2$str[$i]", "$end_final", $conf{c_x}, $conf{c_y}, $alpha_c_0_data2, $delta_c_0_data2, $conf{c_jk});
        $alpha_c_0_data2 = $atmp;
        $delta_c_0_data2 = $dtmp;
    }

# Number of objects
    $ndata1[$i] = count_objects("$data1$str[$i]$end_final");
    $ndata2[$i] = ($ncat == 2 ? count_objects("$data2$str[$i]$end_final") : $ndata1[$i]);

    #$tmp  = `wc -l "$data1$str[$i]$end_final"`;
    #@tmp2 = split(" ", $tmp);
    #$ndata1[$i] = $tmp2[0];
    #die "Error: Wrong number of data points ($ndata1[$i]) for catalogue #$i" unless $ndata1[$i]>1;

    # TODO: write subroutine for extrema

    # Extend of data (for random catalogue)
    if (! -s "$random1$end_final") {
        push(@extr_all1, [extrema("$data1$str[$i]", $conf{c_x}, $conf{c_y})]);
        print "Extrema[z-bin #$i]: $extr_all1[$i][0] $extr_all1[$i][1] $extr_all1[$i][2] $extr_all1[$i][3]\n" unless $quiet;

        if ($i==0) {
            @extr1 = ($extr_all1[0][0], $extr_all1[0][1], $extr_all1[0][2], $extr_all1[0][3]);
        } else {
            $extr1[0] = min($extr1[0], $extr_all1[$i][0]);  # x_min
            $extr1[1] = max($extr1[1], $extr_all1[$i][1]);  # x_max
            $extr1[2] = min($extr1[2], $extr_all1[$i][2]);  # y_min
            $extr1[3] = max($extr1[3], $extr_all1[$i][3]);  # y_max
        }
    }

    if (! -s "$random2$end_final") {
        push(@extr_all2, [extrema("$data2$str[$i]", $conf{c_x}, $conf{c_y})]);
        print "Extrema[z-bin #$i]: $extr_all2[$i][0] $extr_all2[$i][1] $extr_all2[$i][2] $extr_all2[$i][3]\n" unless $quiet;

        if ($i==0) {
            @extr2 = ($extr_all2[0][0], $extr_all2[0][1], $extr_all2[0][2], $extr_all2[0][3]);
        } else {
            $extr2[0] = min($extr2[0], $extr_all2[$i][0]);  # x_min
            $extr2[1] = max($extr2[1], $extr_all2[$i][1]);  # x_max
            $extr2[2] = min($extr2[2], $extr_all2[$i][2]);  # y_min
            $extr2[3] = max($extr2[3], $extr_all2[$i][3]);  # y_max
        }
    }

}

# (Approximate) number of random objects from config entry
$nrand1_max = $conf{nrand};
create_random_cat("$random1", $nrand1_max, $alpha_c_0_data1, $delta_c_0_data1, \@extr1, $conf{c_jk});

# Actual number of random objects
$nrand_outm1 = count_objects("$random1$end_final");

if ($ncat == 2) {
    $nrand2_max = $conf{nrand2} eq "-" ? $nrand1_max : $conf{nrand2};

    create_random_cat("$random2", $nrand2_max, $alpha_c_0_data2, $delta_c_0_data2, \@extr2, $conf{c_jk});
    $nrand_outm2 = count_objects("$random2$end_final");
} else {
    $nrand_outm2 = $nrand_outm1;
}


#############################
### Correlation functions ###
#############################

### Create job lists

# TODO: if $ncat == 2, also run $j < $i

# DD
$index = 0;
foreach $i (0 .. $#str) {
    foreach $j ($i .. $#str) {

        next if defined $options{d} && $i != $j;  # Diagonal only

	$extension[$index] = "D$i" . "D$j";
        $c1[$index]  = "$data1$str[$i]";
        $c2[$index]  = "$data2$str[$j]";
        $nb1[$index] = $conf{nresample};
        $nb2[$index] = $conf{nresample};
        $index ++;
    }
}

# DR (= <D1 R2> for two data catalogues)
foreach $i (0 .. $#str) {
    $extension[$index] = "D$i" . "R";
    $c1[$index]  = "$data1$str[$i]";
    $c2[$index]  = "$random2";
    $nb1[$index] = $conf{nresample};
    # Bootstrap: none for random. Jackknife: same size as data
    $nb2[$index] = $conf{error} eq "bootstrap" ? 0 : $conf{nresample};
    $index ++;
}

# DpR (= <D2 R1> for two data catalogues)
if ($ncat == 2) {
    foreach $i (0 .. $#str) {
        $extension[$index] = "Dp$i" . "R";
        $c1[$index]  = "$data2$str[$i]";
        $c2[$index]  = "$random1";
        $nb1[$index] = $conf{nresample};
	# Bootstrap: none for random. Jackknife: same size as data
        $nb2[$index] = $conf{error} eq "bootstrap" ? 0 : $conf{nresample};
        $index ++;
    }
}

# RR
$extension[$index] = "RR";
$c1[$index]   = "$random1";
$c2[$index]   = "$random2";
if (! defined $options{R}) {
    # New in v1.53
    $nb1[$index] = $conf{nresample};
    $nb2[$index] = $conf{nresample};
} else {
    $nb1[$index] = 0;
    $nb2[$index] = 0;
}
$index ++;
$njob = $index;


### Run correlation jobs

print STDERR "Run 'athena' with $nproc processes in parallel\n" unless $quiet;
init_parallel($nproc, $njob);

# Run threads
foreach $index (0 .. $nproc-1) {
    $tid = run_thread($index, $list_start[$index], $list_stop[$index]);
}
$tid = 0;

# Stop threads
foreach $i (0 .. $nproc-1) {
    $res  = $t[$i]->join();
    die if $res;
}


read_w("w.RR", \@theta, \@w_RR, \@w_RR_resample);
# w_RR_resample unused if jackknife and '-R' options are given.
# In that case: Created identical 'resample' entries from RR for consistency
if (defined $options{R} ) {
    foreach my $b (0 .. $#theta) {
        foreach my $n (0 .. $conf{nresample} - 1) {
            $w_RR_resample[$b][$n] = $w_RR[$b];
        }
    }
}

foreach $i (0 .. $#str) {
# New (v1.3): Factor 2 here instead of in w_LS and w_H
    if ($ncat == 1) {
        $B1[$i] = ($nrand_outm1 - 1)/$ndata1[$i]/2.0;    # <R R>_tot / <Di R>_tot
        $B2[$i] = $B1[$i];                               # <R R>_tot / <Di R>_tot
    } else {
        $B1[$i] = $nrand_outm1 / $ndata1[$i];            # <R1 R2>_tot / <D1i R2>_tot
        $B2[$i] = $nrand_outm2 / $ndata2[$i];            # <R1 R2>_tot / <D2i R1>_tot
        # MKDEBUG: One random catalogue:
        #$B1[$i] = ($nrand_outm1-1) / $ndata1[$i]/2;
        #$B2[$i] = ($nrand_outm2-1) / $ndata2[$i]/2;
    }
}

# Calculate w(theta) and various error contributions
my $nresample = 0;
foreach $i (0 .. $#str) {

    $name_i = "w.D$i" . "R";
    my $tmp = read_w($name_i, \@theta, \@w_DiR, \@w_DiR_resample);
    if ($i == 0) {
        $nresample = $tmp;
    } else {
        die "Number of resample (boots/jack) inconsistent between D0R ($nresample) and D$i" . "R ($tmp)"
            unless $nresample == $tmp;
    }

    foreach $j ($i .. $#str) {

	    next if defined $options{d} && $i != $j;  # Diagonal only

	    $name_ij = "w.D$i" . "D" . "$j";
	    $name_j  = ($ncat == 2 ? "w.Dp$j" . "R" : "w.D$j" . "R");
	    #printf "MKDEBUG $i $j $name_i $name_j\n";
	    $tmp = read_w($name_ij, \@theta, \@w_DiDj, \@w_DiDj_resample);
	    die "Number of resample (boots/jack) inconsistent between D0R ($nresample) and D$i" . "D$j ($tmp)"
	      unless $nresample == $tmp;

	    # Bug in v1.5, produced false result in diagonal mode (option '-d'). Fixed in v1.51
	    # Bug in v1.51, cross-correlations were wrong. Fixed in v1.52
	    # Bug in v1.54, cross-correlations were still wrong! Added $ncat option to read $name_j
	    if ($i != $j || $ncat == 2) {
	      $tmp = read_w($name_j, \@theta, \@w_DjR, \@w_DjR_resample);
	      die "Number of resample (boots/jack) inconsistent between D0R ($nresample) and D$j" . "R ($tmp)"
		unless $nresample == $tmp;
	    } else {
	      @w_DjR = @w_DiR;
	      @w_DjR_resample = @w_DiR_resample;
	    }

	    die "Files $name_ij, $name_j and w.RR have different lengths=#angular bins ($#w_DiDj, $#w_DjR, $#w_RR)+1"
	      if $#w_DiDj!=$#w_DjR || $#w_DjR!=$#w_RR;

	    if ($ncat == 1) {
	        if ($i==$j) {
		        # Auto-correlation. Number of data pairs = D (D-1) / 2
		        $A = $nrand_outm1*($nrand_outm1 - 1)/($ndata1[$i]*($ndata1[$i] - 1));
	        } else {
		        # Cross-correlation. Number of data pairs = Di * Dj
		        $A = $nrand_outm1*($nrand_outm1 - 1)/($ndata1[$i]*$ndata1[$j])/2.0;
	        }
	    } else {
	        # Cross-correlation. Number of data pairs = D1i * D2j, number of random pairs = R1 * R2
	        $A = $nrand_outm1*$nrand_outm2/($ndata1[$i]*$ndata2[$j]);
            # MKDEBUG: One random catalogue:
            #$A = $nrand_outm1*($nrand_outm1-1)/(2*$ndata1[$i]*$ndata2[$j]);
	    }

        # Todo: make sure angular bin values are the same...

	    my $name_base = "w_theta_$i" . "_$j";
	    foreach $estimator ("LS", "Ham") {
	        write_w("$name_base" . "_$estimator.dat", \@w_DiDj, \@w_DiR, \@w_DjR, \@w_RR,
                    $estimator, $A, $B1[$i], $B2[$j], $i, $j, \@theta, \@w_DiDj_resample,
                    \@w_DiR_resample, \@w_DjR_resample, \@w_RR_resample, $nresample, $ncat, 1);

	        write_all_w_resample("$name_base" . "_resample_$estimator.dat", \@w_DiDj_resample,
				     \@w_DiR_resample, \@w_DjR_resample, \@w_RR_resample, $nresample,
				     $estimator, $A, $B1[$i], $B2[$j], $i, $j, \@theta, $ncat)
                if defined $options{r};

	        write_covariance_resample("cov" . "_$i" . "_$j" . "_$estimator.dat", \@w_DiDj_resample,
					  \@w_DiR_resample, \@w_DjR_resample, \@w_RR_resample, $nresample, $estimator,
					  $A, $B1[$i], $B2[$j], $i, $j, $#theta, $ncat, $conf{error})
                unless $conf{error} eq "none";
	    }

        print "Files for $i $j written\n" unless $quiet;

    }

}


# Avoid 'unused' warnings
$nrand_outm1 = $nrand_outm2 = 0.0;
# $#w_RR_resample = 0;



###################
### Subroutines ###
###################

# Returns extrema of coordinates. This might not work for ra/dec!!
sub extrema {
    my ($cat, $cx, $cy) = @_;
    open(CAT, "$cat");

    $xmin = $ymin = 1e30;
    $xmax = $ymax = -$xmin;

    <CAT>;  # first line: nb of objects
    while (<CAT>) {
        @F = split(" ", $_);
        $xmin = $F[$cx] if $xmin>$F[$cx];
        $ymin = $F[$cy] if $ymin>$F[$cy];
        $xmax = $F[$cx] if $xmax<$F[$cx];
        $ymax = $F[$cy] if $ymax<$F[$cy];
    }
    return ($xmin, $xmax, $ymin, $ymax);
}

# Writes the tree code configuration file "config_tree"
sub write_config {
    my ($config_name) = @_;

    my $sformat = defined $conf{c_jk} ? "position_jack_num" : "position";

    open(my $out_fh, ">$config_name");
    print {$out_fh} "GALCAT1         $cat1\n";
    print {$out_fh} "GALCAT2         $cat2\n";
    print {$out_fh} "WCORR           4\n";
    print {$out_fh} "SWCORR_SUBTYPE  nn_2d\n";
    print {$out_fh} "SFORMAT         $sformat\n";
    print {$out_fh} "SCOORD_INPUT    $conf{coord_units}\n";
    print {$out_fh} "SCOORD_OUTPUT   $conf{coord_output}\n";
    print {$out_fh} "THMIN           $conf{thmin}\n";
    print {$out_fh} "THMAX           $conf{thmax}\n";
    print {$out_fh} "NTH             $conf{nth}\n";
    print {$out_fh} "BINTYPE         $conf{bintype}\n";

    my $radec = 1;
    if ($conf{coord_type} eq "xy" || $conf{project} eq "delta_c" || $conf{project} eq "tan") {
        $radec = 0;
    }

    print {$out_fh} "RADEC           $radec\n";
    print {$out_fh} "OATH            $conf{oath}\n";
    print {$out_fh} "SERROR          $conf{error}\n";
    print {$out_fh} "NRESAMPLE       $nresample1 $nresample2\n" unless $conf{error} eq "none";
    close $out_fh;
}

# Runs jobs $start..$stop-1 for thread $i
sub run_thread {
  my ($i, $start, $stop) = @_;

  $t[$i] = async {
    foreach $j ($start .. $stop-1) {
      #print "Thread $i: job $j\n";
      run_job($j);
    }
    return 0;
  };

  return $t[$i]->tid();
}


# Called by run_thread
sub run_job {
  my ($index) = @_;

  #print "MKDBEUG run_job (c1, c2)[$index] = ($c1[$index], $c2[$index]\n" if $index==0;

  run_athena($extension[$index], $c1[$index], $c2[$index], $nb1[$index], $nb2[$index]);
}

# Prepares config file and runs athena
sub run_athena {
  my ($extension, $c1, $c2, $nb1, $nb2) = @_;

  if (! -s "w.$extension") {
    print "### Correlating $extension\n" unless $quiet;

    $cat1        = "$c1$end_final";
    $cat2        = "$c2$end_final";
    $nresample1  = $nb1;
    $nresample2  = $nb2;
    write_config("config_tree.$extension");

    runpr("athena $qflag -c config_tree.$extension --out_w w.$extension --out_log log.$extension");
  } else {
    print STDERR "w.$extension exists, skipping (athena)...\n";
  }
}

# Runs a command and prints it on stderr
# Exits if command does not return 0
sub runpr {
    my ($command) = @_;
    print STDERR "*** Running $athena_bin/$command ***\n";
    $res = `$athena_bin/$command`;

    if ($?!=0) {
       print STDERR "Last command returned error code $?. Stopping "
                    . basename($0) . "\n";
	   exit ($? & 128);
    }

    return $res;
}

# Reads a correlation function file (created by athena)
sub read_w {
    my ($name, $theta, $w, $w_resample) = @_;
    open(W, "$name");
    my $b = 0;
    my $nresample = 0;
    while (<W>) {
        next if /#/;
        @F = split(" ", $_);
        $theta->[$b] = $F[0];
        $w->[$b]     = $F[1];
        foreach $i (2 .. $#F) {
            $w_resample->[$b][$i-2] = $F[$i];
            $nresample = $#F - 1;
        }
        $b++;
    }
    close W;

    return $nresample;
}


### On Gaussian error propagation of the correlation function:
###
### For y = x^n with n=+1,-1:
### Dy^2 = [x^{n-1} Dx]^2 = (y/x Dx)^2
### -> (Dy/y)^2 = (Dx/x)^2
### With Poisson errors: Dx^2 = x
### -> (Dy/y)^2 = 1/|x|
###

# Landy&Szalay 1993 estimator.  No factor 2 in front of the DR-term since DR 
# contains all pairs, DD and RR half of the pairs (no double counting).
sub w_LS {
    my ($DiDj, $DiR, $DjR, $RR, $A, $Bi, $Bj, $i, $j, $verbose) = @_;


    if ($RR == 0) {
        warning("w_LS: No pairs RR found for (i,j)=($i,$j), setting correlation to zero") if defined $verbose;
        return (0, 0);
    }

    my $w = $A*$DiDj/$RR - ($Bi*$DiR + $Bj*$DjR)/$RR + 1;

    my $Bterm = ($Bi**2*$DiR + $Bj**2*$DjR)*2;
    $Bterm = $Bterm / 2.0 if ($i != $j);
    my $dwPsqr = ( $A**2*$DiDj + $Bterm + ($w-1)**2*$RR ) / $RR**2;

    die "Error: Negative error $dwPsqr in w_LS" if $dwPsqr<0;
    my $dwP = sqrt($dwPsqr);

    #printf "MKDEBUG: w_LS(%d, %d) = %.2g * %d/%d - (%.2g * %d + %.2g * %d)/%d + 1\n",
    #  $i, $j, $A, $DiDj, $RR, $Bi, $DiR, $Bj, $DjR, $RR  if $i == $j;

    return ($w, $dwP);
}

# Hamilton 1993 estimator, with factor 4 (2) to account for difference in
# pair numbers between auto (i=j) and cross (i!=j) correlation
sub w_Ham {
    my ($DiDj, $DiR, $DjR, $RR, $i, $j, $ncat, $verbose) = @_;

    if ($DiDj==0) {
        warning("w_Ham: No pairs found (DD) for (i,j)=($i,$j), setting correlation to zero") if defined $verbose;
        return (0, 0);
    }
    if ($DiR==0 || $DjR==0) {
        warning("w_Ham: No pairs (DR) found for (i,j)=($i,$j), setting correlation to zero") if defined $verbose;
        return (0, 0);
    }
    if ($RR == 0) {
        warning("w_Ham: No pairs (RR) found for (i,j)=($i,$j), setting correlation to zero") if defined $verbose;
        return (0,0)
    }

    my $w  = $DiDj*$RR/($DiR*$DjR);

    if ($ncat == 1) {
        $w    *= 2;
        $w    *= 2 if $i == $j;
    }

    $w = $w - 1;

    my $DRterm = 1.0/$DiR + 1.0/$DjR;

    $DRterm *= 2 if $i == $j;

    my $dwPsqr = ($w+1.0)**2 * (1.0/$DiDj + $DRterm + 1.0/$RR);

    die "Error: Negative error $dwPsqr in w_Ham" if $dwPsqr<0;
    my $dwP = sqrt($dwPsqr);

    return ($w, $dwP);
}

# Returns bootstrap or jackknife error
sub rms_resample {
    my ($DiDj, $DiR, $DjR, $RR, $b, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $ncat, $error_type) = @_;

    return (0, 0) if $error_type eq "none";

    my $dw = rms_bootstrap($DiDj, $DiR, $DjR, $RR, $b, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $ncat);

    if ($error_type eq "jackknife") {
	    $dw = $dw / sqrt($nresample) * ($nresample - 1);
    } elsif (! $error_type eq "bootstrap") {
	    die "Invalid error type $error_type";
    }

    return $dw;
}

# Bootstrap variance
sub rms_bootstrap {
    my ($DiDj, $DiR, $DjR, $RR, $b, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $ncat) = @_;

    my $dwB = 0;
    if ($nresample>1) {
	    my $wi = 0;
	    my $wisqr = 0;
	    foreach my $n (0 .. $nresample-1) {
	        if ($estimator eq "LS") {
		        ($x, $dum) = w_LS($DiDj->[$b][$n], $DiR->[$b][$n], $DjR->[$b][$n], $RR->[$b][$n], $A, $Bi, $Bj, $i, $j);
	        } elsif ($estimator eq "Ham") {
		        ($x, $dum) = w_Ham($DiDj->[$b][$n], $DiR->[$b][$n], $DjR->[$b][$n], $RR->[$b][$n], $i, $j, $ncat);
	        }
	        $wi    += $x;
	        $wisqr += $x*$x;
	    }
	    $wi        = $wi/$nresample;
	    my $dwBsqr = ($wisqr - $nresample*$wi*$wi)/($nresample-1);

	    die "Error: Negative error $dwBsqr" if $dwBsqr<0;
	    $dwB = sqrt($dwBsqr);
    } else {
	    $dwB = 0;
    }

    return $dwB;
}

# Returns bootstrap or jackknife covariance element cov_{bc}
sub covariance_resample {
    my ($DiDj, $DiR, $DjR, $RR, $b, $c, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $ncat, $error_type) = @_;

    return 0 if $error_type eq "none";

    my $cov = 0;

    if ($error_type eq "jackknife") {
	    $cov = covariance_jackknife($DiDj, $DiR, $DjR, $RR, $b, $c, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $ncat);
    } elsif ($error_type eq "bootstrap") {
	    $cov = covariance_bootstrap($DiDj, $DiR, $DjR, $RR, $b, $c, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $ncat);
    } else {
	    die "Invalid error type $error_type";
    }

    return $cov;
}

# Jackknife covariance
sub covariance_jackknife {
    my ($DiDj, $DiR, $DjR, $RR, $b, $c, $njack, $estimator, $A, $Bi, $Bj, $i, $j, $ncat) = @_;

    my $cov = 0;
    if ($njack > 1) {
	    my $wib = 0;
	    my $wic = 0;
	    my $wisqr = 0;
	    foreach my $n (0 .. $njack-1) {

	        if ($estimator eq "LS") {
		        ($xb, $dum) = w_LS($DiDj->[$b][$n], $DiR->[$b][$n], $DjR->[$b][$n], $RR->[$b][$n], $A, $Bi, $Bj, $i, $j);
		        ($xc, $dum) = w_LS($DiDj->[$c][$n], $DiR->[$c][$n], $DjR->[$c][$n], $RR->[$c][$n], $A, $Bi, $Bj, $i, $j);
	        } elsif ($estimator eq "Ham") {
		        ($xb, $dum) = w_Ham($DiDj->[$b][$n], $DiR->[$b][$n], $DjR->[$b][$n], $RR->[$b][$n], $i, $j, $ncat);
		        ($xc, $dum) = w_Ham($DiDj->[$c][$n], $DiR->[$c][$n], $DjR->[$c][$n], $RR->[$c][$n], $i, $j, $ncat);
	        }
	        $wib   += $xb;
	        $wic   += $xc;
	        $wisqr += $xb*$xc;
	    }
	    $wib    = $wib/$njack;
	    $wic    = $wic/$njack;
	    $cov = ($wisqr - $njack*$wib*$wic) * ($njack-1) / $njack;

	    die "Error: Negative covariance diagonal $cov ($b, $c)" if $cov<0 && $b == $c;
    } else {
	    $cov = 0;
    }

    return $cov;
}

# Bootstrap covariance
sub covariance_bootstrap {
    my ($DiDj, $DiR, $DjR, $RR, $b, $c, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $ncat) = @_;

    my $cov = 0;
    if ($nresample>1) {
        my $wib = 0;
        my $wic = 0;
        my $wisqr = 0;
        foreach my $n (0 .. $nresample-1) {
            if ($estimator eq "LS") {
                ($xb, $dum) = w_LS($DiDj->[$b][$n], $DiR->[$b][$n], $DjR->[$b][$n], $RR->[$b][$n], $A, $Bi, $Bj, $i, $j);
                ($xc, $dum) = w_LS($DiDj->[$c][$n], $DiR->[$c][$n], $DjR->[$c][$n], $RR->[$c][$n], $A, $Bi, $Bj, $i, $j);
            } elsif ($estimator eq "Ham") {
                ($xb, $dum) = w_Ham($DiDj->[$b][$n], $DiR->[$b][$n], $DjR->[$b][$n], $RR->[$b][$n], $i, $j, $ncat);
                ($xc, $dum) = w_Ham($DiDj->[$c][$n], $DiR->[$c][$n], $DjR->[$c][$n], $RR->[$c][$n], $i, $j, $ncat);
            }
            $wib   += $xb;
            $wic   += $xc;
            $wisqr += $xb*$xc;
        }
        $wib    = $wib/$nresample;
        $wic    = $wic/$nresample;
        $cov = ($wisqr - $nresample*$wib*$wic)/($nresample-1);

        # Note: If $nresample-1 is replaced by $nresample, the resulting covariance matches the one from
        # columns2cov.pl + covariance.pl

        die "Error: Negative covariance diagonal $cov ($b, $c)" if $cov<0 && $b == $c;
    } else {
        $cov = 0;
    }

    return $cov;
}

# Write angular correlation function with Poisson and bootstrap error bars to a file
sub write_w {
    my ($name, $w_DiDj, $w_DiR, $w_DjR, $w_RR, $estimator, $A, $Bi, $Bj, $i, $j, $theta,
            $w_DiDj_resample, $w_DiR_resample, $w_DjR_resample, $w_RR_resample, $nresample, $ncat, $verbose) = @_;

    open(my $w_fh, ">$name");
    print {$w_fh} "# theta[$conf{coord_output}]              w    err_resample     err_Poisson. Estimator = $estimator, resample = $conf{error}\n";

    my $w = 0;
    my $dw_Poi = 0;
    foreach my $b (0 .. $#theta) {

	if ($estimator eq "LS") {
	    ($w, $dw_Poi) = w_LS($w_DiDj->[$b], $w_DiR->[$b], $w_DjR->[$b], $w_RR->[$b], $A, $Bi, $Bj, $i, $j, $verbose);
	} elsif ($estimator eq "Ham") {
	    ($w, $dw_Poi) = w_Ham($w_DiDj->[$b], $w_DiR->[$b], $w_DjR->[$b], $w_RR->[$b], $i, $j, $ncat, $verbose);
	}
	my $dw_B = rms_resample($w_DiDj_resample, $w_DiR_resample, $w_DjR_resample, $w_RR_resample, $b, $nresample,
				$estimator, $A, $Bi, $Bj, $i, $j, $ncat, $conf{error});

	printf {$w_fh} "%12.8f   % 15.6g %15.6g %15.6g\n", $theta->[$b], $w, $dw_B, $dw_Poi;

    }

    close $w_fh;
}

# Writes the angular correlation functions from all resamples to file "$name". Called if option '-r' is set.
sub write_all_w_resample {
    my ($name, $w_DiDj_resample, $w_DiR_resample, $w_DjR_resample, $w_RR_resample, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $theta, $ncat) = @_;

    open(my $w_fh, ">$name");
    print {$w_fh} "# theta[$conf{coord_output}]       w_{i=0..Nresample}, estimator = $estimator\n";

    foreach my $b (0 .. $#theta) {

	printf {$w_fh} "%12.8f  ", $theta->[$b];
	foreach my $n (0 .. $nresample - 1) {
	    if ($estimator eq "LS") {
	        ($w_resample, $dum) = w_LS($w_DiDj_resample->[$b][$n], $w_DiR_resample->[$b][$n], $w_DjR_resample->[$b][$n],
					   $w_RR_resample->[$b][$n], $A, $Bi, $Bj, $i, $j);
	    } elsif ($estimator eq "Ham") {
		($w_resample, $dum) = w_Ham($w_DiDj_resample->[$b][$n], $w_DiR_resample->[$b][$n], $w_DjR_resample->[$b][$n],
					$w_RR_resample->[$b][$n], $i, $j, $ncat);
	    }
	    printf {$w_fh} " %15.6g", $w_resample;
	}
	print {$w_fh} "\n";
	$dum = 0;
    }

    close $w_fh;
}

# Writes the resampled (bootstrap or jackknife) covariance to a file
sub write_covariance_resample {
    my ($name, $w_DiDj_resample, $w_DiR_resample, $w_DjR_resample, $w_RR_resample, $nresample, $estimator,
	$A, $Bi, $Bj, $i, $j, $Ntheta, $ncat, $error_type) = @_;

    open(my $cov_fh, ">$name");

    my @cov = ();
    foreach my $b (0 .. $Ntheta) {
	foreach my $c ($b .. $Ntheta) {
	    $cov[$b][$c] = covariance_resample($w_DiDj_resample, $w_DiR_resample, $w_DjR_resample, $w_RR_resample,
					       $b, $c, $nresample, $estimator, $A, $Bi, $Bj, $i, $j, $ncat, $error_type);
	}
    }

    # Symmetrize
    foreach my $b (0 .. $Ntheta) {
	foreach my $c (0 .. $b-1) {
	    $cov[$b][$c] = $cov[$c][$b];
	}
    }

    foreach my $b (0 .. $Ntheta) {
	foreach my $c (0 .. $Ntheta) {
	    printf {$cov_fh} "%15.12g ", $cov[$b][$c];
	}
	print {$cov_fh} "\n";
    }

    close $cov_fh;
}


# Puts columns in right order and projects coordinates if necessary
sub create_final_cat {

    my ($name, $end_final, $cx, $cy, $alpha_c_in, $delta_c_in, $cjk) = @_;

    my $name_final = "$name$end_final";

    $jk_flag = defined $cjk ? "-jk $cjk" : "";

    if ($conf{coord_type} eq "radec" && $conf{project} ne "none") {

        # radec -> xy (projection)

        if (!defined $delta_c_in) {
	        # Calculates center of data points for projection
	        $tmp = runpr("center_gal.pl -ra $cx -dec $cy -nohead $name");
	        @tmp2 = split(" ", $tmp);
	        $alpha_c = $tmp2[0];
	        $delta_c = $tmp2[1];
	        die "Coordinate center of catalogue $name could not be calculated" unless defined $delta_c;
	        print "(alpha_c, delta_c) = ($alpha_c, $delta_c)\n" unless $quiet;
	    } else {
	        $alpha_c = $alpha_c_in;
	        $delta_c = $delta_c_in;
	        print "Using input (alpha_c, delta_c) = ($alpha_c, $delta_c)\n" unless $quiet;
	    }

	    # Project to Cartesian coordinates
	    if (! -s $name_final) {
	        runpr("cat2gal.pl $qflag -coord_input $conf{coord_units} -coord_output $conf{coord_units} -x $cx -y $cy $jk_flag -alphac $alpha_c -deltac $delta_c $name tan > $name_final");
	    } else {
	        print "$name_final exists, skipping (cat2gal)...\n";
	    }

    } else {

	    if ($conf{coord_type} eq "radec") {

	        # radec -> radec (no Transformation)
	        if (! -s $name_final) {
                runpr("cat2gal.pl $qflag -coord_input $conf{coord_units} -coord_output $conf{coord_units} -x $cx -y $cy $jk_flag $name no > $name_final");
	        } else {
                print "$name_final exists, skipping (cat2gal)...\n";
	        }

	    } elsif ($conf{coord_type} eq "xy") {

	        #  xy -> xy (no Transformation)
	        if (! -s $name_final) {
                runpr("cat2gal.pl $qflag -coord_input $conf{coord_units} -coord_output $conf{coord_units} -x $cx -y $cy $jk_flag $name no > $name_final");
	        } else {
                print "$name_final exists, skipping (cat2gal)...\n";
	        }
	    }

    }

    return ($alpha_c, $delta_c);
}

# Create random catalogue (call venice)
sub create_random_cat {
    my ($random, $nrand_max, $alpha_c_0, $delta_c_0, $extr, $c_jk) = @_;

    if (! -s "$random$end_final") {
	    if (! -s $random) {
	        $coord_in = $conf{coord_type} eq "radec" ? "spher" : "cart";
	        if ($conf{mask} eq "-") {
		        # No mask given (e.g. simulated data)
		        runpr("venice -r -xmin $extr->[0] -xmax $extr->[1] -ymin $extr->[2] -ymax $extr->[3] -o $random -npart $nrand_max -coord $coord_in");
	        } else {
		        # MKDEBUG: Make sure, patched version with clock/time gsl init is compiled
		        runpr("venice -r -m $conf{mask} -o $random -npart $nrand_max -coord $coord_in");
	        }
	    } else {
	        print "$random exists, skipping (venice)...\n";
	    }
	    create_final_cat("$random", "$end_final", 0, 1, $alpha_c_0, $delta_c_0, $c_jk);
    } else {
	    print "$random$end_final exists, skipping (create_final_cat)...\n";
    }
}

# Runs 'wc' to count objects in the file '$cat'
sub count_objects {
    my ($cat) = @_;

    my $tmp  = `wc -l "$cat"`;
    my @tmp2 = split(" ", $tmp);
    my $n    = $tmp2[0];
    return $n;
}

# Reads the values from a configuration file into the hash %conf
sub read_config {
    my ($filename) = @_;

    open(CONFIG, $filename);

    print "Reading config file $filename\n";
    while (<CONFIG>) {
        next if /^\s*#/;
        @F = split(" ", $_);
        next if $#F==-1;

        # Scalar entries
        foreach $key (@conf_entry) {
            $conf{"$key"} = $F[2] if $F[0] =~ /\b$key\b/ && $F[1] eq "=";
        }

        # Vector entries
        if ($F[0] =~ /specifiers/) {
            for $z (2 ..$#F) {
                # Stop if comment starts
                last if $F[$z] =~ /#/;
                push(@specifiers, $F[$z]);
            }
        }

        # Optional entries
        foreach $key (@conf_entry_opt) {
            $conf{"$key"} = $F[2] if $F[0] =~ /\b$key\b/ && $F[1] eq "=";
        }

        # Undefine entries that are -1
        foreach $key (@conf_entry_opt) {
            $conf{"$key"} = undef if defined $conf{"$key"} and $conf{"$key"} == -1;
        }
    }
    close CONFIG;
}


# Exits if invalid config value is found
sub check_config {

    #local $SIG{__DIE__}; 

    foreach $key (@conf_entry) {
        if (! defined $conf{"$key"}) {
            die "*** Entry '$key' not found in config_woftheta ***\n";
        }
    }

    die "Key 'coord_type' in config_woftheta has to be 'radec' or 'xy'." 
        if !($conf{coord_type} eq "radec" || $conf{coord_type} eq "xy");

    die "Key 'project' in config_woftheta has to be 'none', 'cosdelta_c' or 'tan'" 
	if !($conf{project} eq "none" || $conf{project} eq "cosdelta_c" || $conf{project} eq "tan");

    die "Combination coord=\"xy\" and project!=\"none\" does not make sense"
        if $conf{coord_type} eq "xy" && $conf{project} ne "none";

    die "Error has to be 'none', 'bootstrap' or 'jackknife'" unless
        $conf{error} eq "none" or $conf{error} eq "bootstrap" or $conf{error} eq "jackknife";

    die "nresample cannot be zero unless error is 'none'" if
        ! ($conf{error} eq "none") and $conf{nresample} == 0;

    die "nresample has to be zero if error is 'none'" if
        $conf{nresample} != 0 and $conf{error} eq "none";

    die "c_jk (Jackknife sample number column) is defined but nresample is zero" if
        defined $conf{c_jk} and $conf{nresample} == 0;

}

# Prints the contents of the config file, read by read_config()
sub out_config {
    print STDERR "*** Config info: ***\n";
    foreach $key (keys %conf) {
	    print STDERR "$key -> $conf{$key}\n";
    }
    print STDERR "specifiers -> ";
    for $z (0 .. $#specifiers) {
	    print STDERR "$specifiers[$z] ";
    }
    print "\n";
    print STDERR "*** End config info ***\n";
}

# Prints a warning
sub warning {
    my ($message) = @_;
    print STDERR "Warning: $message\n" unless $quiet;
}

# Cleans up from a previous run (interactive)
sub clean_up {

  my ($all) = @_;

  print STDERR "Clean-up from previous run:\n";

  # RR
  if (-e "w.RR") {
    print "Delete 'w.RR' (if 'y', 'athena' will be re-run to obtain RR correlation)? [n] ";
    if (! defined $all) {
      chomp($ans = <STDIN>);
    } else {
      print "y\n"; $ans = "y";
    }
    unlink "w.RR" if $ans eq "y";
  }
  
  # DR
  $#list = -1;
  foreach $i (0 .. $#str) {
    $end_i = "D$i" . "R";
    unshift @list, "w.$end_i" if -e "w.$end_i";
  }
  if ($ncat == 2) {
    foreach $i (0 .. $#str) {
      $end_i = "Dp$i" . "R";
      unshift @list, "w.$end_i" if -e "w.$end_i";
    }
  }
  if ($#list!=-1) {
    print "Delete 'w.D*R' (if 'y', 'athena' will be re-run to obtain DiR correlations)? [n] ";
    if (! defined $all) {
      chomp($ans = <STDIN>);
    } else {
      print "y\n"; $ans = "y";
    }
    unlink @list if $ans eq "y";
  }

  # DD
  $#list = -1;
  foreach $i (0 .. $#str) {
    foreach $j ($i .. $#str) {
      $end_ij = "D$i" . "D$j";
      unshift @list, "w.$end_ij" if -e "w.$end_ij";
    }
  }
  if ($#list!=-1) {
    print "Delete 'w.D*D*' (if 'y', 'athena' will be re-run to obtain DiDj correlations)? [n] ";
    if (! defined $all) {
      chomp($ans = <STDIN>);
    } else {
      print "y\n"; $ans = "y";
    }
    unlink @list if $ans eq "y";
  }

  # Formatted data
  $#list = -1;
  foreach $i (0 .. $#str) {
    unshift @list, "$data1$str[$i]$end_final" if -e  "$data1$str[$i]$end_final";
  }
  if ($#list!=-1) {
    print "Delete '$data1*$end_final' (if 'y', 'cat2gal' will be re-run to format data catalogue)? [n] ";
    if (! defined $all) {
      chomp($ans = <STDIN>);
    } else {
      print "y\n"; $ans = "y";
    }
    unlink @list if $ans eq "y";
  }

  # Formatted data2
  if ($ncat == 2) {
      $#list = -1;
      foreach $i (0 .. $#str) {
	  unshift @list, "$data2$str[$i]$end_final" if -e  "$data2$str[$i]$end_final";
      }
      if ($#list!=-1) {
	  print "Delete '$data2*$end_final' (if 'y', 'cat2gal' will be re-run to format data catalogue)? [n] ";
	  if (! defined $all) {
	      chomp($ans = <STDIN>);
	  } else {
	      print "y\n"; $ans = "y";
	  }
	  unlink @list if $ans eq "y";
      }
  }

  # Formatted random
  if (-e "$random1$end_final") {
    print "Delete '$random1$end_final' (if 'y', 'cat2gal' will be re-run to format random catalogue)? [n] ";
    if (! defined $all) {
      chomp($ans = <STDIN>);
    } else {
      print "y\n"; $ans = "y";
    }
    unlink "$random1$end_final" if $ans eq "y";
  }

  # Formatted random2
  if ($ncat == 2) {
      if (-e "$random2$end_final") {
	  print "Delete '$random2$end_final' (if 'y', 'cat2gal' will be re-run to format random catalogue)? [n] ";
	  if (! defined $all) {
	      chomp($ans = <STDIN>);
	  } else {
	      print "y\n"; $ans = "y";
	  }
	  unlink "$random2$end_final" if $ans eq "y";
      }
  }

  # Random
  if (-e "$random1") {
    print "Delete '$random1' (if 'y', 'venice' will be re-run to obtain random catalogue)? [n] ";
    if (! defined $all) {
      chomp($ans = <STDIN>);
    } else {
      print "y\n"; $ans = "y";
    }
    unlink "$random1" if $ans eq "y";
  }

  # Random2
  if ($ncat == 2) {
      if (-e "$random2") {
          print "Delete '$random2' (if 'y', 'venice' will be re-run to obtain random catalogue)? [n] ";
          if (! defined $all) {
              chomp($ans = <STDIN>);
          } else {
              print "y\n"; $ans = "y";
          }
          unlink "$random2" if $ans eq "y";
      }
  }

}

# Initialises job lists for the distribution on multiple CPUs
sub init_parallel {
    my ($ncpu, $njob) = @_;

    use integer;
    my $naverage = $njob / $ncpu;  # Zero if $ncpu>$njob
        my $nextra   = $njob % $ncpu;
    no integer;

    print "#cpus = $ncpu, #jobs = $njob\n" unless $quiet;
    print "#average = $naverage, #extra = $nextra\n" unless $quiet;

    my $start = 0;
    my $stop = 0;
    foreach $i (0 .. $ncpu-1) {

# Create job list for cpu #i
        $stop += $naverage;
        $stop++ if $i<$nextra;

        if ($start>=$stop) {
            $list_start[$i] = $list_stop[$i] = -1;
        } else {
            $list_start[$i] = $start;
            $list_stop[$i]  = $stop;
        }

        $start = $stop;
    }
}

# Prints help string
sub usage {
    my ($ex) = @_;

    print STDERR "Usage: woftheta_xcorr.pl [OPTIONS]\n";
    print STDERR "OPTIONS:\n";
    print STDERR "   -c CONFIG    Config file CONFIG (default: 'config_woftheta')\n";
    print STDERR "   -n           Do not query and do not clean up from a previous run\n";
    print STDERR "   -a           Delete all relevant files from a previous run without query\n";
    print STDERR "   -p NPROC     Run 'athena' with NPROC processes in parallel (default 1)\n";
    print STDERR "   -d           'Diagonal' only (no cross-correlation between different z)\n";
    print STDERR "   -r           Ouput not only resample errors and covariance but all resample\n";
    print STDERR "                 correlation functions (files 'w_theta_i_j_resample_{LS,Ham}.dat'\n";
    print STDERR "   -R           For Jackknife errors, do not resample RR pairs\n";
    print STDERR "   -q           Quiet(er) mode\n";
    print STDERR "   -h           This message\n";

    exit $ex if defined $ex;
}

