#!/usr/bin/env perl
###############################################################################
#
#    barcodeByCU.pl
#    
#    Input is a "squished" file of contigs, IE. all non-coding regions are
#    removed and both strands are orientated in the forward direction. The
#    length of all input sequences MUST be a multiple of 3. This script will
#    calculate the codon usage for each contig in the input file.
#
#    Copyright (C) Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;
 
#CPAN modules
use Bio::SeqIO;
use Data::Dumper;
use Bio::Tools::CodonTable;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands 
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
my $global_options = checkParams();
if(!exists $global_options->{'silent'}) { printAtStart(); }

######################################################################
# CODE HERE
######################################################################
# override defaults
my $global_cut_off_len = overrideDefault(102, 'cutoff');
my $global_protein_code =  overrideDefault(11,'protein');

# make the kmer array
my %global_prot_hash = ();              # reverse lookup of proteins
my @global_mer_array = makeMers();      # store all the 3-mers

# keep count of where we're at
my $global_update_amount = 1000;

my $global_valid_seqs = 0;
my $global_parsed_seqs = 0;
my $global_complete_seqs = 0;
my $global_rejected_seqs = 0;

# open the output file
my $global_out_fh = openWrite($global_options->{'out'});
printOutHeaderRow($global_out_fh);

# one thread for parsing through the seqio object
my $seqio = Bio::SeqIO->new(-file => $global_options->{'in'}, '-format' => 'Fasta');

while(my $seq = $seqio->next_seq) {
    # record seq length
    my $seq_len = $seq->length();
    next if ($seq_len < $global_cut_off_len);
  
    my @raw_barcode = cutMers4Barcode($seq->seq);
    foreach my $protein (keys %global_prot_hash)
    {
        my $total = 0;
        foreach my $index (@{$global_prot_hash{$protein}}) 
        {
            $total += $raw_barcode[$index];
        }
        foreach my $index (@{$global_prot_hash{$protein}}) 
        {
            if(0 != $total)
            {
                $raw_barcode[$index] = sprintf("%0.4f", ($raw_barcode[$index]/$total));
            }
        }
    }
    print $global_out_fh join "\t", @raw_barcode, "\n";    
}

# close the files
close($global_out_fh);

if(!exists $global_options->{'silent'})
{
print<<EOF
    Processed: $global_complete_seqs sequences
    Rejected: $global_rejected_seqs sequences
----------------------------------------------------------------
EOF
}

######################################################################
# CUSTOM SUBS
######################################################################
sub printOutHeaderRow
{
    #-----
    # print the header row for the csv files
    #
    my ($fh) = @_;
    print $fh "\"SequenceID\"";
    foreach my $kmer (@global_mer_array)
    {
        print $fh ",\"$kmer\"";
    }
    print $fh "\n";
}

sub cutMers4Barcode {
    #-----
    # cut each window in kmers and do barcodes
    #
    my ($sequence) = @_;
    my $sequence_length = length $sequence;
    my $sub_start = 0;
    
    # get a fresh map to do counting in
    my $mer_map_ref = merArray2Map();
    
    # cut into kmers and add to the map
    while($sub_start + 3 <= $sequence_length)
    {
        # look out for non ACGT chars!
        my $this_mer = substr $sequence, $sub_start, 3;
        print "$this_mer\n"; 
        if(exists $$mer_map_ref{$this_mer})
        {
            $$mer_map_ref{$this_mer}++;
        }
        $sub_start += 3;
    }
    
    my @raw_barcode = ();
    foreach my $this_mer (@global_mer_array)
    {
        push @raw_barcode, ($$mer_map_ref{$this_mer});
    }
    return @raw_barcode;
}

sub merArray2Map {
    #-----
    # return an map made from the @global_mer_array
    #
    my %return_map = ();
    foreach my $this_mer (@global_mer_array)
    {
        $return_map{$this_mer} = 0;
    }
    return \%return_map;
}

sub makeMers {
    #-----
    # Make all 3-mers and work out their groupings
    #
    my @alphabet = ('A', 'C', 'G', 'T');
    my @output_array = ();
    my $index = 0;
    my $codon_table  = Bio::Tools::CodonTable -> new ( -id => $global_protein_code );
    foreach my $i (0..3) {
        foreach my $j (0..3) {
            foreach my $k (0..3) {
                my $mer = $alphabet[$i].$alphabet[$j].$alphabet[$k];
                my $protein  = $codon_table->translate($mer);
                if(!exists $global_prot_hash{$protein})
                {
                    # first time seen - get an array to store all the relative mers
                    my @tmp = ();
                    push @tmp, $index;
                    $global_prot_hash{$protein} = \@tmp;
                } 
                else
                {
                    push @{$global_prot_hash{$protein}}, $index;
                }
                
                # push the trimer onto the main array
                push @output_array, $mer;
                
                # update the index so we can store the rev-lookup right
                $index++;
            }
        }
    }
    return @output_array;
}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS
  
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help+", "in|i:s", "out|o:s", "cutoff|c:i", "silent+", "protein|p:i",);
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    if(!exists $options{'in'} ) { printParamError ("We need an input file to process!"); }
    if(!exists $options{'out'} ) { printParamError ("We need an output file to process!"); }
    #if(!exists $options{''} ) { printParamError (""); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #  
    my ($error) = @_;  
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

######################################################################
# FILE IO
sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          }, 
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    # 
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map { $_ . " " . $params->{$_}} keys %{$params};
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well    
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : "  . $! . "\n";
        }
    }
}

######################################################################
# MISC

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

  barcodeByCU.pl
  
=head1 COPYRIGHT

   copyright (C) Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

  Input is a "squished" file of contigs, IE. all non-coding regions are removed
  and both strands are orientated in the forward direction. The length of all
  input sequences MUST be a multiple of 3. This script will calculate the codon usage
  for each contig in the input file.
  
=head1 SYNOPSIS

    barcodeByCU.pl -in|i FILE -out|o FILE

    -i -in FILE                  The SQUISHED file to work on
    -o -out FILE                 The file to write barcodes to
    [-c -cutoff INT]             Reject all sequences shorter than this [default: 102]
    [-protein -p CODON_CODE]     Protein translation table [default: 11] --> see below 
    [-silent]                    Output nothing extra to the screen
    [-help]                      Displays basic usage information

    CODON_CODE
     
    Specify a number from the following list (Uses: Bio::Tools::CodonTable)
      
    1 Standard
    2 Vertebrate Mitochondrial
    3 Yeast Mitochondrial
    4 Mold, Protozoan,_and_CoelenterateMitochondrial_and_Mycoplasma/Spiroplasma
    5 Invertebrate Mitochondrial
    6 Ciliate, Dasycladacean_and_Hexamita_Nuclear
    9 Echinoderm Mitochondrial
    10 Euplotid Nuclear
    11 Bacterial
    12 Alternative Yeast_Nuclear
    13 Ascidian Mitochondrial
    14 Flatworm Mitochondrial
    15 Blepharisma Nuclear
    16 Chlorophycean Mitochondrial
    21 Trematode Mitochondrial
    22 Scenedesmus obliquus_Mitochondrial
    23 Thraustochytrium Mitochondrial
         
=cut
