use strict;
use Getopt::Std;
use Bio::SeqIO;

sub help_mess{
	print "\t", "-i <File>", "\t", "Input fasta file", "\n";
    print "\t", "-h [Bool]", "\t", "Help. Show this and exits.", "\n";
}

## ARGS #############
my %options;
my $arguments = 'i:8h';
getopts($arguments, \%options);
my $_input = $options{'i'};
my $_debug = $options{'8'};
my $_help = $options{'h'};
%options = ();
if($_help){help_mess();exit;}

my $io = Bio::SeqIO->new(-file => $_input, -format => 'fasta'); 
while(my $seq_obj = $io->next_seq){
    my $seq = $seq_obj -> seq;
    my $id = $seq_obj -> display_id;
    while($seq =~ /([nN]+)/g){
        print "$id:$-[0]-", $+[0]-1, "\t", length($1), "\n";
    }
}
