#!/usr/bin/env perl

# gibbs_sampler_motifs.pl v0.1.0
# Author: MH Seabolt
# Last updated: 2022-03-22

# SYNOPSIS:
# This program just contains example driver code to instantiate a Gibbs sampler object and 
# pass it some DNA sequences to search for potential motifs of interest.  
# Extend this program and teach it new tricks however you like :)


##################################################################################
# The MIT License
#
# Copyright (c) 2022 Matthew H. Seabolt
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
my $version;
my $help;

sub usage {
	my $usage = "gibbs_sampler_motifs.pl v0.1.0\n
	PURPOSE: This program just contains example driver code to instantiate a Gibbs sampler object and 
		 pass it some DNA sequences to search for potential motifs of interest. 
			
	USAGE:	gibbs_sampler_motifs.pl -i <sequences.fasta> -o <motif.profile.tab> -k <3> -l <7>
	
	INPUT/OUTPUT:
	-i | --input 		FILE; input sequences in FASTA format (sequences are expected to be the same length, but dont necessarily have to be aligned)
	-o | --output		STR; output filename
	
	SAMPLING PARAMETERS:
	-k | --kiter		INT; Number of iterations to check for convergence ( Default: 3 )
	-l | --len 		INT; Length of the motif you want to sample for ( Default: 7 -- minimum value 5 )
	-n | --nrepet		INT; Number of replicates to sample from a random starting point ( Default: 100 )
	
	OPTIONAL EXTRAS:
	-v | --version	 	Print version number and exit.
	-h | --help		Print this help message and exit.
	
	\n";
	print $usage;
}

GetOptions(	'input|i=s'  => \$input,
			'out|o=s'    => \$output,
			'k|kiter=i'  => \$k,
			'len|l=i'    => \$motif_len,
			'nrepet|n=i' => \$nsamples,
			'version|v'  => \$version,
			'help|h'     => \$help,
) or die usage();

# Print the version number or the help message and exit if -v or -h is activated
if ( $version ) 	{ die "gibbs_sampler_motifs.pl v0.1.0\n"; 		}
if ( $help    )     { die usage();											}

# Sanity check the other user-defined parameters
$k = ( $k && int($k) >= 2 )? $k : 3;
$motif_len = ( $motif_len && int($motif_len) >= 5 )? $motif_len : 7;
$nsamples = ( $nsamples && int($nsamples) >= 1 )? $nsamples : 100;

# Parse the given FASTA sequences into an array
my $DNA = read_fasta($input);

# Instantiate a new Gibbs object
# Note that the array of sequences should be passed as a reference
my $Sampler = Gibbs->new( seqs => $DNA, k => $k, motif_len => $motif_len );

# Use our Gibbs object to run the MCMC sampling procedure.
# The returned value is a tuple with the motif sequence and its score.
my $best_motif = $Sampler->sample( $nsamples );
my ( $motif, $score ) = @{ $best_motif };
print STDOUT "\n\nBest motif found: $motif       Score: $score  \n";

# Print the converged profile out if we have an output file specified, otherwise just print to STDOUT
if ( $output ne "--" )	{	$Sampler->print_profile( $output, ">" );		}
else 					{	$Sampler->print_profile();						}		# Prints to STDOUT


# The arguments here are given to generalize the computation using any size data, but will default to the specifics of the $Sampler object if params are omitted.
# Here we are determining the minimum length the sequences should be if we want to reasonably rule out the likelihood that a random, unrelated 6-mer 
# occurs in a total of 400 sequences by random chance (using a standard ATCG alhpabet for DNA).  Our threshold for "reasonable likelihood" is 99% surety (aka p <= 0.01)
# Order of arguments here is ( probability threshold, number of sequences, motif length, alphabet size )
my $min_sequence_length = $Sampler->minimum_suggested_sequence_length( 0.01, 400, 6, 4 );
print STDOUT "Minimum suggested sequence length: $min_sequence_length\n";

# Determine the probability of a random 6-mer occuring in all 400 DNA sequences if we know all the other parameters.
# This is essentially the inverse of the above problem -- here we know the details of the sequences, but want to determine the likelihood of "noisy" results
my $reasonable_likelihood = $Sampler->nontarget_motif_probability(400, 6, 18314, 4);
print STDOUT "Likelihood of a random 6-mer appearing in all 400 DNA sequences of length 18314 bp: $reasonable_likelihood\n";

exit;

############# SUBROUTINES ###################

# This is just an example subroutine to basically parse a FASTA file
#     and return an array of sequences, which is what we need to create our Gibbs object.
#     You can create any sort of parsing routines that you like, as long as the end 
#     result is an array of sequences in the same upper or lower case (preferably upper)
sub read_fasta 		{
	my ($input) = @_;
	
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
	
	return $DNA;
}
