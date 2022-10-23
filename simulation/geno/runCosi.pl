#!/usr/bin/perl -w
use strict;

my $path_to_cosi = $ARGV[0];
my $locus_number = $ARGV[1];
my $pop1_index = $ARGV[2];
my $pop2_index = $ARGV[3];
my $recParam_file = $ARGV[4];
my $param_file = $ARGV[5];
my $out_prefix = $ARGV[6];

my $rd = 1;

my $ttl_rnds = 100;

while ($rd <= $ttl_rnds)
{
	print "\nRound $rd:\n";
	my $length;
	
	open(IN, $param_file) || die "Could not open params file.\n";
	while (<IN>) {
	    if (m/^length/) {
		$length = $';
		chomp($length);
		$length =~ s/\s+//g;
	    }
	}
	
	print "$length\n\n";
	print "$recParam_file\n";
	my $resp = `$path_to_cosi/cosi_package/recosim $recParam_file $length`;
	print $resp;
	$resp = `$path_to_cosi/cosi_package/coalescent -p $param_file -o $out_prefix`;
	print $resp;
	
	my $outpos_pop1 = $out_prefix.".pos-".$pop1_index;
	my $outhap_pop1 = $out_prefix.".hap-".$pop1_index;
	my $outpos_pop2 = $out_prefix.".pos-".$pop2_index;
        my $outhap_pop2 = $out_prefix.".hap-".$pop2_index;

	my $newpos_pop1 = $out_prefix."_round".$rd.".pos-".$pop1_index;
	my $newhap_pop1 = $out_prefix."_round".$rd.".hap-".$pop1_index;
	my $newpos_pop2 = $out_prefix."_round".$rd.".pos-".$pop2_index;
        my $newhap_pop2 = $out_prefix."_round".$rd.".hap-".$pop2_index;
	
	#if(!(-e $newpos_pop1)) {system "mv $outpos $newpos";}
	#else {die "$newpos exits\n";}
	
	#if(!(-e $newhap)) {system "mv $outhap $newhap";}
	#else {die "$newhap exits\n";}
	
	$rd ++;

}

print "\n";
