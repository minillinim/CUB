#!/usr/bin/env perl
###############################################################################
#
#    cub.pl
#    
#    Barcode a set of contigs by codon usage
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
printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################
my $global_gff_file = $global_options->{'in'}.".gff";
my $global_orfs_fasta_file = $global_gff_file.".fasta";
my $global_reject_length = overrideDefault(50,'length');

# first call the orfs on the contigs file
checkAndRunCommand("callOrfs.pl", [{
                                  -in => $global_options->{'in'},
                                  -out => $global_gff_file
                                  }], DIE_ON_FAILURE);

# make the gff3 file into a multiple fasta
checkAndRunCommand("gff2fasta.pl", [{
                                   -gff => $global_gff_file,
                                   -fasta => $global_options->{'in'},
                                   -out => $global_orfs_fasta_file,
                                   -w => "0",
                                   -l => $global_reject_length
                                   }], DIE_ON_FAILURE);

# update the orfs file ready for squishing
my $pre_squished_file = $global_orfs_fasta_file;
if(exists $global_options->{'keep'})
{
    # the user wishes to keep the orfs
    # we need to squish elsewhere
    $pre_squished_file = $global_orfs_fasta_file.".presquish";
    checkAndRunCommand("sed", [{
                              -e => "\"s/_[^_]*_[^_]*_.$//\""
                              },
                              "$global_orfs_fasta_file > $pre_squished_file"
                              ], DIE_ON_FAILURE);  
}
else
{
    checkAndRunCommand("sed", [
                              -i,
                              {
                              -e => "\"s/_[^_]*_[^_]*_.$//\""
                              },
                              $pre_squished_file
                              ], DIE_ON_FAILURE);  
}

my $squished_file = $pre_squished_file.".squished";
# squish it all up       
checkAndRunCommand("squishOrfs.pl", [{
                                    -in => $pre_squished_file,
                                    -out => $squished_file
                                    }], DIE_ON_FAILURE);
   
# clean up the ORFs file (if need be)   
checkAndRunCommand("rm", [{}, $pre_squished_file], WARN_ON_FAILURE);

# now barcode!             
checkAndRunCommand("barcodeByCU.pl", [{
                                     -i => $squished_file,
                                     -o => $global_options->{'out'},
                                     -t => $global_options->{'threads'},
                                     -c => $global_reject_length
                                     }], DIE_ON_FAILURE); 

# remove the squished file
checkAndRunCommand("rm", [{}, $squished_file], WARN_ON_FAILURE); 

######################################################################
# CUSTOM SUBS
######################################################################

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "in|i:s", "out|o:s", "glimmer|g+", "keep|k+", "length|l:i");
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
    if(!exists $options{'in'} ) { printParamError ("It is theoretically impossible to barcode contigs if you do not supply a contigs file..."); }
    if(!exists $options{'out'} ) { printParamError ("Please tell me where to write the contigs to"); }
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
    my $param_str = join " ", map {formatParams($_)} @{$params};
    
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to 
    # checkAndRunCommand
    #
    my $ref = shift;
    
    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
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
            carp "**WARNING: $0 : " . $! . "\n";
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

    cub.pl

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

   Barcode a set of contigs by codon usage

=head1 SYNOPSIS

    cub.pl -in|i CONTIGS_FILE -out|o BARCODE_FILE

      -in -i CONTIGS_FILE          Contigs file to barcode
      -out -o BARCODE_FILE         Output file
      [-length -l LENGTH]          Reject any orfs shorter than this length [default: 50]
      [-keep -k]                   Keep the ORFs file
      [-glimmer -g]                Use glimmer to call orfs [default: use prodigal]
      [-help -h]                   Displays basic usage information
         
=cut