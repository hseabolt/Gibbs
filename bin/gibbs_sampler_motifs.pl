#!/usr/bin/perl

# gibbs_sampler_motifs.pl v0.1.0
# Author: MH Seabolt
# Last updated: 1-13-2020

# SYNOPSIS:
# This program just contains example driver code to instantiate a Gibbs sampler object and 
# pass it some DNA sequences to search for potential motifs of interest.  Extend
# this program and teach it new tricks however you like :)

##################################################################################
# The MIT License
#
# Copyright (c) 2021 Matthew H. Seabolt
#
# Permission is hereby granted, free of charge, 
# to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to 
# deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom 
# the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice 
# shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##################################################################################


use strict;
use warnings;
use Getopt::Long qw(GetOptions);


# Load our Gibbs package
use lib '/scicomp/home/ngr8/Biolinux/scripts_for_release/Gibbs/bin';
use Gibbs;

# Required input parameters
my $input = "--";
my $output = "--";
my $k;
my $motif_len;
my $nsamples;

sub usage {
	my $usage = "gibbs_sampler_motifs.pl v0.1.0\n
	PURPOSE: This program just contains example driver code to instantiate a Gibbs sampler object and 
			pass it some DNA sequences to search for potential motifs of interest. Extend
			this program and teach it new tricks however you like :)
			
	USAGE:	gibbs_sampler_motifs.pl -i sequences.fasta -o motif.profile.tab -k 3 -l 7
	
	INPUT/OUTPUT:
	-i	input sequences in FASTA format (sequences are expected to be the same length, but dont necessarily have to be aligned)
	-o 	output file name, including extensions
	
	SAMPLING PARAMETERS:
	-k	INT; Number of iterations to check for convergence ( Default: 3 )
	-l	INT; Length of the motif you want to sample for ( Default: 7 -- minimum value 5 )
	-n 	INT; Number of replicates to sample from a random starting point ( Default: 100 )
	\n";
	print $usage;
}

GetOptions(	'input|i=s' => \$input,
			'out|o=s' => \$output,
			'k=i' => \$k,
			'len|l=i' => \$motif_len,
			'nrepet|n=i' => \$nsamples,
) or die usage();

# Parameter Setups
$k = ( $k && int($k) >= 2 )? $k : 3;
$motif_len = ( $motif_len && int($motif_len) >= 5 )? $motif_len : 5;
$nsamples = ( $nsamples && int($nsamples) >= 1 )? $nsamples : 100;

# Read the sequences into an array
$/ = ">";
my $fh = *STDIN;
my $succin = open(DATA, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *DATA if ( $succin ); 
	my @fastas = <$fh>;
	my $trash = shift @fastas;	# Get rid of the first element, which will be a lone ">" symbol
close DATA if ( $succin );
$/ = "\n";

my $DNA = [ ];
foreach my $record ( @fastas )	{
	my ($header, @seq) = split "\n", $record;
	my $sequence = join '', @seq;
	$sequence =~ s/>//g;						# Remove any lingering ">" symbols
	push @{$DNA}, uc $sequence;
}

# Instantiate a new Gibbs object
my $Sampler = Gibbs->new( DNA => $DNA, k => $k, motif_len => $motif_len );
my $best_motif = $Sampler->sample_random_starting_positions( $nsamples );
print STDERR "\n\nBest motif found: $best_motif->[0]       Score: $best_motif->[1]\n";

# Print the converged profile out if we have an output file specified, otherwise just print to STDOUT
if ( $output ne "--" )	{	$Sampler->print_profile( $Sampler->{_profiles}->[-1], ">", "$output" );		}
else 					{	$Sampler->print_profile( $Sampler->{_profiles}->[-1] );						}		# Prints to STDOUT


exit;

