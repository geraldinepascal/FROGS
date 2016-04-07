#!/usr/bin/perl

#    FASTX-toolkit - FASTA/FASTQ preprocessing tools. fastx_barcode_splitter.pl
#    Copyright (C) 2009  A. Gordon (gordon@cshl.edu)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
use strict;
use warnings;
use IO::Handle;
use IO::Zlib;
use PerlIO::gzip;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Pod::Usage;

=pod

=head1 NAME
    
    splitbc.pl

=head1 SYNOPSIS
    
	splitbc.pl read_file_r1 [read_file_r2] --bcfile FILE --prefix-r1 r1.%.fq [--prefix-r2 r2.%.fq] --bol|--eol] 
         [--mismatches N] [--exact] [--partial N] [--trim|--trim2] [--no_adapt] [--rad STR --radTAG STR] [--help] [--quiet] [--debug]

=head1 DESCRIPTION
    
    splitbc.pl classify fasta/fastq single or paired end reads in function of barcode forward or reverse in the first or both reads.
    Also, whether reads are RAD sequence, it checks that barcode is followed by RAD tag.

=head1 OPTIONS

=over
=item B<--bcfile FILE>

	Barcodes file name. (splitbc.pl --doc for details)
	
=item B<--prefix-r1 PREFIX>

	File prefix. The sample name will be placed where the % is placed.
	
=item B<--prefix-r2 PREFIX>

	File prefix. The sample name will be placed where the % is placed.

=item B<--bol>

	Try to match barcodes at the BEGINNING of sequences.
	(What biologists would call the 5' end, and programmers
	 would call index 0.)
	 
=item B<--eol>

	Try to match barcodes at the END of sequences.
	(What biologists would call the 3' end, and programmers
	would call the end of the string.). Give barcode as they appear in the sequence file, so reverse complement if needed.
	NOTE: one of --bol, --eol must be specified, but not both.

=item B<--mismatches N>

	Max. number of mismatches allowed. default is 1.

=item B<--exact>

	Same as '--mismatches 0'. If both --exact and --mismatches are specified, '--exact' takes precedence.

=item B <--partial N>
	
	Allow partial overlap of barcodes. (splitbc.pl --doc for details)
	(Default is not partial matching)
	
item B<--trim>
	
	Should the barecode be trimmed.

=item B<--trim2>

	Shoud the read 2 be trimmed to have the same length as the read1
	NOTE: one of --trim, --trim2 must be specified, but not both.

=item B<--no_adapt>
	
	there is no adaptator (T or A ) between barcode and the sequence

=item B<--rad STR>

	sequence are RADSeq, barcode is immediately followed by rad sequence <STR>.
	NOTE: --rad implies --no_adapt

=item B<--radTAG STR>

	sequence retain at the begining of the sequence after cliping.
	NOTE: if define rad, radTAG must be defined to, and conversely
	Barcode will be splited if they are folowed by radTAG, otherwise they will be recorded in the unmatched file.
	
=item B<--TAG_mismatch N>

	Max. number of mismatches allowed in the radTAG sequence. default is 0.

=item B<--quiet>

	Don't print counts and summary at the end of the run.
	(Default is to print.)

=item B<--debug>

	Print lots of useless debug information to STDERR.

=item B<--help>

	This helpful help screen.

=back

=head1 AUTHORS

	Jérôme Mariette
    modified by Maria Bernard

=head1 VERSION

    1.1

=head1 DATE

    modified 21/05/2013
    
=head1 KEYWORDS

    Barcode RAD

=head1 EXAMPLE

	(Assuming 's_2_100.txt' is a FASTQ file, 'mybarcodes.txt' is the barcodes file):
   $0 cat s_2_100.txt | $0 --bcfile mybarcodes.txt --bol --mismatches 2 \\
   	--prefix-r1 /tmp/bla_%.txt"
        
=head1 CHANGELOG

	1.1: 18/10/2013 Maria
		construct list of best barcode to deal with equal match. In this case reads will be written in "ambiguous" files
           
=cut

##
## This program splits a FASTQ/FASTA file into several smaller files,
## Based on barcode matching.
##
## run with "--help" for usage information
##
## Assaf Gordon <gordon@cshl.edu> , 11sep2008

# Forward declarations
sub load_barcode_file ($);
sub parse_command_line ;
sub match_sequences ;
sub mismatch_count($$) ;
sub print_results;
sub open_and_detect_input_format;
sub read_record;
sub write_record($);
sub usage();

# Global flags and arguments, 
# Set by command line argumens
my $barcode_file ;
my $barcodes_at_eol = 0 ;
my $barcodes_at_bol = 0 ;
my $exact_match = 0 ;
my $allow_partial_overlap = 0;
my $allowed_mismatches = 1;
my $newfile_prefix_r1 = '';
my $newfile_prefix_r2 = '';
my $quiet = 0 ;
my $debug = 0 ;
my $trim = 0 ;
my $trim2 = 0 ;
my $adapt = 1 ;
my $no_adapt = 0 ;
my $rad;
my $radTAG;
my $TAG_mm=0;
my $fastq_format = 1;

# Global variables 
# Populated by 'parse_command_line'
my $radTAG_length=0;
my $check_rad=0;

# Global variables 
# Populated by 'create_output_files'
my %filenames;
my %files;
my %counts = ( 'unmatched' => 0 );
my $barcodes_length;
my @barcodes;
my $input_file_io;
my $input_file_io_r2;
my $paired = 0;
my $write_r2 = 0;

# The Four lines per record in FASTQ format.
# (when using FASTA format, only the first two are used)
my $seq_name;
my $seq_bases;
my $seq_name2;
my $seq_qualities;
my $seq_len;
my $seq_name_r2;
my $seq_bases_r2;
my $seq_name2_r2;
my $seq_qualities_r2;
my $seq_len_r2;

#
# Start of Program
#
parse_command_line ;

load_barcode_file ( $barcode_file ) ;

open_and_detect_input_format;

match_sequences ;

print_results unless $quiet;

#
# End of program
########################################################################


sub parse_command_line {
	my $help;
	my $doc; 
	my $version;

	pod2usage(1) if (scalar @ARGV==0);

	my $result = GetOptions ( "bcfile=s" => \$barcode_file,
				  "eol"  => \$barcodes_at_eol,
				  "bol"  => \$barcodes_at_bol,
				  "exact" => \$exact_match,
				  "prefix-r1=s" => \$newfile_prefix_r1,
				  "prefix-r2=s" => \$newfile_prefix_r2,
				  "quiet" => \$quiet, 
				  "partial=i" => \$allow_partial_overlap,
				  "trim" => \$trim,
				  "trim2" => \$trim2,
				  "debug" => \$debug,
				  "mismatches=i" => \$allowed_mismatches,
				  "no_adapt"=> \$no_adapt,
				  "rad=s"=> \$rad,
				  "radTAG=s"=> \$radTAG,
				  "TAG_mismatch=i"=> \$TAG_mm,
				  "help" => \$help,
				  "doc" => \$doc,
				  "version" => \$version
				  ) ;
	
	pod2usage(1) if ($help);
	doc() if ($doc);
	version() if ($version);
	
	die "Error: barcode file not specified (use '--bcfile [FILENAME]')\n" unless defined $barcode_file;
	die "Error: prefix path/filename not specified (use '--prefix-r1 [PATH]%.fq')\n" unless defined $newfile_prefix_r1;
	# If a read2 file is provided 
	if ($#ARGV == 1) {
		die "Error: prefix path/filename not specified (use '--prefix-r2 [PATH]%.fq')\n" unless defined $newfile_prefix_r2;
	}
	
	if ($barcodes_at_bol == $barcodes_at_eol) {
		die "Error: can't specify both --eol & --bol\n" if $barcodes_at_eol;
		die "Error: must specify either --eol or --bol\n" ;
	}

	die "Error: invalid for value partial matches (valid values are 0 or greater)\n" if $allow_partial_overlap<0;

	$allowed_mismatches = 0 if $exact_match;

	die "Error: invalid value for mismatches (valid values are 0 or more)\n" if ($allowed_mismatches<0);

	die "Error: partial overlap value ($allow_partial_overlap) bigger than " . 
		"max. allowed mismatches ($allowed_mismatches)\n" if ($allow_partial_overlap > $allowed_mismatches);

	$trim = 1 if $trim2 ;
	
	if ((defined $rad && !defined $radTAG) || (!defined $rad && defined $radTAG)){
		die "Error: you must defined --rad AND --radTAG option\n";
	}
	if ( $no_adapt){
		$adapt=0;
	}
	if (defined $rad ){
		chomp $rad;
		chomp $radTAG;
		die "Error: $rad does not ended by $radTAG\n" if ($rad !~ /$radTAG$/g);
		$adapt=0;
		$check_rad=1;
		#~ warn $rad, " ",$radTAG;
		$radTAG_length=length($radTAG);
	}
	exit unless $result;
}



#
# Read the barcode file
#
sub load_barcode_file ($) {
	my $filename = shift or croak "Missing barcode file name";

	open BCFILE,"<$filename" or die "Error: failed to open barcode file ($filename)\n";
	while (<BCFILE>) {
		next if m/^#/;
		chomp $_;
		my @arr = split(' ',$_);
		my $ident = $arr[0];
		my $barcode = uc($arr[1]);

		# Sanity checks on the barcodes
		die "Error: bad data at barcode file ($filename) line $.\n" unless defined $barcode;
		die "Error: bad barcode value ($barcode) at barcode file ($filename) line $.\n"
			unless $barcode =~ m/^[AGCT]+$/;

		die "Error: bad identifier value ($ident) at barcode file ($filename) line $. (must be alphanumeric)\n" 
			unless $ident =~ m/^\w+$/;

		die "Error: badcode($ident, $barcode) is shorter or equal to maximum number of " .
		    "mismatches ($allowed_mismatches). This makes no sense. Specify fewer  mismatches.\n" 
		    	if length($barcode)<=$allowed_mismatches;

		$barcodes_length = length($barcode) unless defined $barcodes_length;
		die "Error: found barcodes in different lengths. this feature is not supported yet.\n" 
			unless $barcodes_length == length($barcode);

	 	push @barcodes, [$ident, $barcode];

		if ($allow_partial_overlap>0) {
			foreach my $i (1 .. $allow_partial_overlap) {
				substr $barcode, ($barcodes_at_bol)?0:-1, 1, '';
	 			push @barcodes, [$ident, $barcode];
			}
		}
	}
	close BCFILE;

	if ($debug) {
		print STDERR "barcode\tsequence\n";
		foreach my $barcoderef (@barcodes) {
			my ($ident, $seq) = @{$barcoderef};
			print STDERR $ident,"\t", $seq ,"\t",length($seq),"\n";
		}
	}
}

# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub create_output_files {
	my %barcodes = map { $_->[0] => 1 } @barcodes; #generate a uniq list of barcode identifiers;
	$barcodes{'unmatched'} = 1 ;
	$barcodes{'ambiguous'} = 1 ;
	
	foreach my $ident (keys %barcodes) {
		
		my $new_filename = $newfile_prefix_r1;
		$new_filename =~ s/%/$ident/g;
		$filenames{$ident} = $new_filename;
		open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n"; 
		$files{$ident} = $file ;
		
		if (defined $rad && $ident ne "unmatched" && $ident ne "ambiguous"){
			my $new_filename = $newfile_prefix_r1;
			$new_filename =~ s/%/${ident}_2rad/g;
			$filenames{$ident."_2rad"} = $new_filename;
			open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n"; 
			$files{$ident."_2rad"} = $file ;	
		}
		
		if ($paired == 1) {
			my $new_filename_r2 = $newfile_prefix_r2;
			$new_filename_r2 =~ s/%/$ident/g;
			$filenames{$ident."r2"} = $new_filename_r2;
			open my $file, ">$new_filename_r2" or die "Error: failed to create output file ($new_filename)\n"; 
			$files{$ident."r2"} = $file ;
			
			if (defined $rad && $ident ne "unmatched" && $ident ne "ambiguous"){
				my $new_filename_r2 = $newfile_prefix_r2;
				$new_filename_r2 =~ s/%/${ident}_2rad/g;
				$filenames{$ident."_2radr2"} = $new_filename_r2;
				open my $file, ">$new_filename_r2" or die "Error: failed to create output file ($new_filename)\n"; 
				$files{$ident."_2radr2"} = $file ;	
			}
		}
	}
}

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub match_sequences {

	my %barcodes = map { $_->[0] => 1 } @barcodes; #generate a uniq list of barcode identifiers;
	$barcodes{'unmatched'} = 1 ;
	$barcodes{'ambiguous'} = 1 ;
	
	#reset counters
	foreach my $ident ( keys %barcodes ) {
		$counts{$ident} = 0;
	}

	create_output_files;

	# Read file FASTQ file
	# split accotding to barcodes
	do {
		chomp $seq_bases;
		#~ print STDERR "sequence ".length($seq_bases)." $seq_bases: \n" if $debug;
		if ($fastq_format){
			chomp $seq_qualities;
		}
		
		if ($paired == 1) {
			chomp $seq_bases_r2;
			if ($fastq_format){
				chomp $seq_qualities_r2;
			}
		}
		my $best_barcode_mismatches_count = $barcodes_length;
		my $best_barcode_ident = undef;
		my $sequence_fragment;
		my $sequence_fragment_radTAG;
		my $start = 0;
		my $cut_len = $seq_len;
		my $start2 = 0;
		my $cut_len2 = $seq_len_r2;
#		warn $start , " " , $cut_len , " " , $start2 , " " , $cut_len2;
		#Try all barcodes, find the one with the lowest mismatch count
		if ($barcodes_at_bol) {
			$sequence_fragment = substr $seq_bases, 0, $barcodes_length;
			if (defined $rad)
			{	
				$sequence_fragment_radTAG = substr $seq_bases, $barcodes_length, $radTAG_length;
			}
			if ($trim){
				$start=$barcodes_length+$adapt;
				$cut_len=$seq_len-$barcodes_length-$adapt;
				if ($trim2) {
					$start2 = 0 ;
					$cut_len2=$cut_len;
				}
			}
		} elsif($barcodes_at_eol) {
			if ($paired != 1){
				$sequence_fragment = substr $seq_bases, - $barcodes_length ;
				if (defined $rad){	
					$sequence_fragment_radTAG = substr $seq_bases, - ($barcodes_length+$adapt+$radTAG_length), $radTAG_length;
				}
				if ($trim){
					$start=0 ;
					$cut_len=$seq_len-$barcodes_length-$adapt;
				}
			} else {
				$sequence_fragment = substr $seq_bases_r2, - $barcodes_length;
				if (defined $rad){	
					$sequence_fragment_radTAG = substr $seq_bases_r2, - ($barcodes_length+$adapt+$radTAG_length), $radTAG_length;
				}
				if ($trim){
					$start2=0 ;
					$cut_len2=$seq_len_r2-$barcodes_length-$adapt;
					if ($trim2) {
						$start = 0 ;
						$cut_len=$cut_len2 ;
					}
				}				
			}	
			#~ warn "sequence fragment : $sequence_fragment\n" if $debug;
			#~ $sequence_fragment = reverse_complement $sequence_fragment;
			#~ warn "sequence fragment revcomp: $sequence_fragment \n" if $debug;
		}
#		warn $start," ", $cut_len, " ",$start2, " ",$cut_len2;

		foreach my $barcoderef (@barcodes) {
			my ($ident, $barcode) = @{$barcoderef};

			# Get DNA fragment (in the length of the barcodes)
			# The barcode will be tested only against this fragment
			# (no point in testing the barcode against the whole sequence)
			my $mm=$barcodes_length;
			
			# check validity of RAD.
			if (defined $rad){	
				my $mm_rad = mismatch_count($sequence_fragment_radTAG, $radTAG) ; 
				if ($mm_rad<= $TAG_mm){ $mm = mismatch_count($sequence_fragment, $barcode) ; }
			}else {
				$mm = mismatch_count($sequence_fragment, $barcode) ; 
			}

			# if this is a partial match, add the non-overlap as a mismatch
			# (partial barcodes are shorter than the length of the original barcodes)
			$mm += ($barcodes_length - length($barcode)); 

			if ( $mm < $best_barcode_mismatches_count ) {
				$best_barcode_mismatches_count = $mm ;
				$best_barcode_ident = $ident ;
			}elsif($mm == $best_barcode_mismatches_count){
				$best_barcode_ident = "ambiguous"
			}
		}
		
		if ($best_barcode_mismatches_count<=$allowed_mismatches && $best_barcode_ident ne 'ambiguous'){ 
			$seq_bases = substr $seq_bases, $start, $cut_len;
			if ($fastq_format) {
				$seq_qualities = substr $seq_qualities, $start, $cut_len;
			}
			if ($paired == 1 ) {
				$seq_bases_r2 = substr $seq_bases_r2, $start2, $cut_len2;
				if ($fastq_format) {
					$seq_qualities_r2 = substr $seq_qualities_r2, $start2, $cut_len2;
				}
			}
		}
		
		if ( (!defined $best_barcode_ident) || ($best_barcode_mismatches_count>$allowed_mismatches)){ 
			$best_barcode_ident = 'unmatched'; 
		}

		#get the file associated with the matched barcode.
		#(note: there's also a file associated with 'unmatched' barcode)
		
		if (defined $rad && $best_barcode_ident ne 'unmatched' && $best_barcode_ident ne 'ambiguous'){
			my $nb_rad =()= ($seq_bases =~ /$rad/g);
			#~ warn $seq_bases," ",$nb_rad;
			$best_barcode_ident.='_2rad' if ($nb_rad > 0);
			#~ warn $best_barcode_ident;
		}
		
		$counts{$best_barcode_ident}++;
		my $file = $files{$best_barcode_ident};
		#~ warn "$seq_name , $best_barcode_ident \n";
		#~ warn $file ; 
		$write_r2 = 0;
		write_record($file);
		if ($paired == 1) {
			#~ warn $best_barcode_ident."r2";
			my $file = $files{$best_barcode_ident."r2"};
			$write_r2 = 1;
			#~ warn $best_barcode_ident."r2";
			write_record($file);			
		}
		
	} while ( read_record );
}

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { 
	length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) 
}



sub print_results
{
	print "Barcode\tCount\tLocation\n";
	my $total = 0 ;
	foreach my $ident (sort keys %counts) {
		print $ident, "\t", $counts{$ident};
		if ($paired) {
			print "(*2)";
		}		
		print "\t",$filenames{$ident};
		if ($paired) {
			print " & ".$filenames{$ident."r2"};
		}
		print "\n";
		$total += $counts{$ident};	
	}
	print "total\t",$total;
	if ($paired) {
		print "(*2)";
	}
	print "\n";
}


sub read_record
{
	$seq_name = $input_file_io->getline();

	return undef unless defined $seq_name; # End of file?

	$seq_bases = $input_file_io->getline();
	die "Error: bad input file, expecting line with sequences\n" unless defined $seq_bases;
	$seq_len = length($seq_bases)-1;
	# If using FASTQ format, read two more lines
	if ($fastq_format) {
		$seq_name2  = $input_file_io->getline();
		die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2;

		$seq_qualities = $input_file_io->getline();
		die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities;
	}
	
	if ($paired == 1) {
		$seq_name_r2 = $input_file_io_r2->getline();

		return undef unless defined $seq_name_r2; # End of file?
	
		$seq_bases_r2 = $input_file_io_r2->getline();
		die "Error: bad input file, expecting line with sequences\n" unless defined $seq_bases_r2;
		$seq_len_r2 = length($seq_bases_r2)-1;	
		# If using FASTQ format, read two more lines
		if ($fastq_format) {
			$seq_name2_r2  = $input_file_io_r2->getline();
			die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2_r2;
	
			$seq_qualities_r2 = $input_file_io_r2->getline();
			die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities_r2;
		}
	}
	
	return 1;
}

sub write_record($)
{
	if ($write_r2 == 1) {
		my $file = shift;
		croak "Bad file handle" unless defined $file;
	
		print $file $seq_name_r2;
		print $file $seq_bases_r2,"\n";
	
		#if using FASTQ format, write two more lines
		if ($fastq_format) {
			print $file $seq_name2_r2;
			print $file $seq_qualities_r2,"\n";
		}
	} else {
		my $file = shift;
		croak "Bad file handle $seq_name" unless defined $file;
	
		print $file $seq_name;
		print $file $seq_bases,"\n";
	
		#if using FASTQ format, write two more lines
		if ($fastq_format) {
			print $file $seq_name2;
			print $file $seq_qualities,"\n";
		}
	}
}

sub open_and_detect_input_format
{
#	$input_file_io  = new IO::Handle;
#	die "Failed to open STDIN " unless $input_file_io->fdopen(fileno(STDIN),"r");

	if ($ARGV[0] =~ /.gz$/){
		die "Failed to open '$ARGV[0]' " unless open ($input_file_io, "<:gzip", $ARGV[0] ) ;
		#~ $input_file_io = IO::Zlib->new();
		#~ die "Failed to open '$ARGV[0]' " unless $input_file_io->open($ARGV[0],"r");
	}else{
		$input_file_io= IO::Handle->new();
		die "Failed to open '$ARGV[0]' " unless open($input_file_io , $ARGV[0]);
	}
	
	
	# If a read2 file is provided 
	if ($#ARGV == 1) {
		if ($ARGV[1] =~ /.gz$/){
			die "Failed to open '$ARGV[0]' " unless open ($input_file_io_r2 , "<:gzip", $ARGV[1] ) ;
			#~ $input_file_io_r2 = IO::Zlib->new();
			#~ die "Failed to open '$ARGV[0]' " unless $input_file_io_r2->open($ARGV[1],"r");
		}else{
			$input_file_io_r2 = IO::Handle->new();
			die "Failed to open '$ARGV[1]' " unless open($input_file_io_r2 , $ARGV[1]);
		}
		$paired = 1;
	}

	# Get the first characeter, and push it back
	my $first_char = $input_file_io->getc();
	$input_file_io->ungetc(ord $first_char);
	$seq_name = $input_file_io->getline();
	$seq_bases = $input_file_io->getline();
	$seq_len = length($seq_bases)-1;
	if ($paired == 1) {
		$seq_name_r2 = $input_file_io_r2->getline();
		$seq_bases_r2 = $input_file_io_r2->getline();
		$seq_len_r2 = length($seq_bases_r2)-1;		
	}
	if ($first_char eq '>') {
		# FASTA format
		$fastq_format = 0 ;
		print STDERR "Detected FASTA format\n" if $debug;
	} elsif ($first_char eq '@') {
		# FASTQ format
		$fastq_format = 1;
		$seq_name2  = $input_file_io->getline();
		die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2;
		$seq_qualities = $input_file_io->getline();
		die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities;
		if ($paired == 1) {
			$seq_name2_r2  = $input_file_io_r2->getline();
			die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2_r2;
			$seq_qualities_r2 = $input_file_io_r2->getline();
			die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities_r2;			
		}
		print STDERR "Detected FASTQ format\n" if $debug;
	} else {
		die "Error: unknown file format. First character = '$first_char' (expecting > or \@)\n";
	}
#	warn "len1: ", $seq_len, "; len2: ", $seq_len_r2, " ;barcode len: ", $barcodes_length, "; no adapt: ", $adapt;
}

sub doc()
{
print<<EOF;
Barcode Splitter, by Assaf Gordon (gordon\@cshl.edu), 11sep2008

This program reads FASTA/FASTQ file and splits it into several smaller files,
Based on barcode matching.
FASTA/FASTQ data is read from STDIN (format is auto-detected.)
Output files will be writen to disk.
Summary will be printed to STDOUT.

usage: $0 read_file_r1 [read_file_r2] --bcfile FILE --prefix-r1 r1.%.fq [--prefix-r2 r2.%.fq] --bol|--eol] 
         [--mismatches N] [--exact] [--partial N] [--trim|--trim2] [--no_adapt] [--rad STR --radTAG STR] [--help] [--quiet] [--debug]

Arguments:

--bcfile FILE	- Barcodes file name. (see explanation below.)
--prefix-r1 PREFIX	- File prefix. The sample name will be placed where the % is placed.
--prefix-r2 PREFIX	- File prefix. The sample name will be placed where the % is placed.
--bol		- Try to match barcodes at the BEGINNING of sequences.
		  (What biologists would call the 5' end, and programmers
		  would call index 0.)
--eol		- Try to match barcodes at the END of sequences.
		  (What biologists would call the 3' end, and programmers
		  would call the end of the string.). Give barcode as they appear in the sequence file, so reverse complement if needed.
		  NOTE: one of --bol, --eol must be specified, but not both.
--mismatches N	- Max. number of mismatches allowed. default is 1.
--exact		- Same as '--mismatches 0'. If both --exact and --mismatches 
		  are specified, '--exact' takes precedence.
--partial N	- Allow partial overlap of barcodes. (see explanation below.)
		  (Default is not partial matching)
--trim		- Should the barecode be trimmed.
--trim2		- Shoud the read 2 be trimmed to have the same length as the read1
		  NOTE: one of --trim, --trim2 must be specified, but not both.
--no_adapt	- there is no adaptator (T or A ) between barcode and the sequence
--rad STR	- sequence are RADSeq, barcode is immediately followed by rad sequence <STR>.
		  NOTE: --rad implies --no_adapt
--radTAG STR	- sequence retain at the begining of the sequence after cliping.
		  NOTE: if define rad, radTAG must be defined to, and conversely
		  Barcode will be splited if they are folowed by radTAG, otherwise they will be recorded in the unmatched file.
--TAG_mismatch N	- Max. number of mismatches allowed in the radTAG sequence. default is 0.	  
--quiet		- Don't print counts and summary at the end of the run.
		  (Default is to print.)
--debug		- Print lots of useless debug information to STDERR.
--help		- This helpful help screen.

Example (Assuming 's_2_100.txt' is a FASTQ file, 'mybarcodes.txt' is 
the barcodes file):

   \$ cat s_2_100.txt | $0 --bcfile mybarcodes.txt --bol --mismatches 2 \\
   	--prefix-r1 /tmp/bla_%.txt"

Barcode file format
-------------------
Barcode files are simple text files. Each line should contain an identifier 
(descriptive name for the barcode), and the barcode itself (A/C/G/T), 
separated by a TAB character. Example:

    #This line is a comment (starts with a 'number' sign)
    BC1 GATCT
    BC2 ATCGT
    BC3 GTGAT
    BC4 TGTCT

For each barcode, a new FASTQ file will be created (with the barcode's 
identifier as part of the file name). Sequences matching the barcode 
will be stored in the appropriate file.

Running the above example (assuming "mybarcodes.txt" contains the above 
barcodes), will create the following files:
	/tmp/bla_BC1.txt
	/tmp/bla_BC2.txt
	/tmp/bla_BC3.txt
	/tmp/bla_BC4.txt
	/tmp/bla_unmatched.txt
The 'unmatched' file will contain all sequences that didn't match any barcode.

Barcode matching
----------------

** Without partial matching:

Count mismatches between the FASTA/Q sequences and the barcodes.
The barcode which matched with the lowest mismatches count (providing the
count is small or equal to '--mismatches N') 'gets' the sequences.

Example (using the above barcodes):
Input Sequence:
    GATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG

Matching with '--bol --mismatches 1':
   GATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
   GATCT (1 mismatch, BC1)
   ATCGT (4 mismatches, BC2)
   GTGAT (3 mismatches, BC3)
   TGTCT (3 mismatches, BC4)

This sequence will be classified as 'BC1' (it has the lowest mismatch count).
If '--exact' or '--mismatches 0' were specified, this sequence would be 
classified as 'unmatched' (because, although BC1 had the lowest mismatch count,
it is above the maximum allowed mismatches).

Matching with '--eol' (end of line) does the same, but from the other side
of the sequence.

** With partial matching (very similar to indels):

Same as above, with the following addition: barcodes are also checked for
partial overlap (number of allowed non-overlapping bases is '--partial N').

Example:
Input sequence is ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
(Same as above, but note the missing 'G' at the beginning.)

Matching (without partial overlapping) against BC1 yields 4 mismatches:
   ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
   GATCT (4 mismatches)

Partial overlapping would also try the following match:
   -ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
   GATCT (1 mismatch)

Note: scoring counts a missing base as a mismatch, so the final
mismatch count is 2 (1 'real' mismatch, 1 'missing base' mismatch).
If running with '--mismatches 2' (meaning allowing upto 2 mismatches) - this 
seqeunce will be classified as BC1.

EOF

exit 1;
}

sub version()
{
	print STDOUT "1.0";
	exit 0 ;
}