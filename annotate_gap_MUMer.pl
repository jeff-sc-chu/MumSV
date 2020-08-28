use strict;
use Getopt::Std;

sub help_mess{
	print "\t", "-i <File>", "\t", "Mummer diff file", "\n";
    print "\t", "-j [File]", "\t", "Mummer coord file", "\n";
    print "\t", "-g [File]", "\t", "Reference genome Gap info file", "\n";
	print "\t", "-G [File]", "\t", "Query gneome Gap info file", "\n";
    print "\t", "-c [File]", "\t", "chrom name correspondence for Ref", "\n";
	print "\t", "-C	[File]", "\t", "chrom name correspondence for Query", "\n";
    print "\t", "-h [Bool]", "\t", "Help. Show this and exits.", "\n";
}

## ARGS #############
my %options;
my $arguments = 'i:o:j:J:c:C:g:G:8h';
getopts($arguments, \%options);
my $_input = $options{'i'};
my $_output = $options{'o'};
my $_debug = $options{'8'};
my $_help = $options{'h'};
my $_coord_file = $options{'j'};
my $_rGapFile = $options{'g'};
my $_qGapFile = $options{'G'};
my $_refNames = $options{'c'};
my $_queryNames = $options{'C'};
%options = ();
if($_help){help_mess();exit;}
if(!$_output){
	open OUT, ">&STDOUT";
}
else{
	open OUT, "+>$_output";
}	

if(!$_coord_file){die "Error: Required coord file missing.\nUse -h to see help information\n";}
if(!$_rGapFile || !$_qGapFile){die "Error: Required Gap file missing.\nUse -h for help\n";}

my %nameTable;
populateNames();
AnnotateRGap();
AnnotateQGap();


sub AnnotateRGap{
	my %gap;
	open IN, $_rGapFile;
	while(<IN>){
		chomp;
		my ($coord, $size) = split /\t/;
		my ($c, $s, $e) = $coord =~ /(\S+):(\d+)-(\d+)/;
		$gap{$c}{$s} = $e;
	}
	close IN;
	
	foreach my $file($_input, $_coord_file){
		open IN, $file;
		my($prefix, $suffix) = $file =~ /(.*)\.(\w+)/;
		open OUT, "+>$prefix.rGapAnno.$suffix";
		my @starts;
		my $previous_chr;
		while(my $l = <IN>){
			chomp $l;
			my @info = split /\t/, $l;
			my ($c, $s, $e);
			if($file eq $_input){
				$c = $info[0];
				($s, $e) = $info[2] < $info[3] ? ($info[2], $info[3]) : ($info[3], $info[2]);
			}
			elsif($file eq $_coord_file){
				$c = $info[13];
				($s, $e) = $info[0] < $info[1] ? ($info[0], $info[1]) : ($info[0], $info[1]);
			}
			my $gapString;
			if($c ne $previous_chr){
				@starts = sort {$a<=>$b} keys %{$gap{$c}};
				$previous_chr = $c;
			}
			foreach my $gapS(@starts){
				last if($gapS > $e);
				if(isOverlap($s, $e, $gapS, $gap{$c}{$gapS})){
					$gapString .= "RN:$gapS-$gap{$c}{$gapS};";
				}
			}
			print OUT $gapString ? "$l\t$gapString\n" : "$l\t.;\n";
		}
		close IN;
	}	
}

sub AnnotateQGap{
	my %gap;
	open IN, $_qGapFile;
	while(<IN>){
		chomp;
		my ($coord, $size) = split /\t/;
		my ($c, $s, $e) = $coord =~ /(\S+):(\d+)-(\d+)/;
		$c = $nameTable{$c};
		$gap{$c}{$s} = $e;
	}
	close IN;
	
	my($prefix, $suffix) = $_input =~ /(.*)\.(\w+)/;
	open IN, "$prefix.rGapAnno.$suffix";
	open OUT, "+>$prefix.gapAnno.$suffix";
	my @starts;
	my $previousChr;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my ($c, $s, $e, $cAlt);
		my ($invC, $invS, $invE);
		my ($previousBlock, $postBlock) = getBlocks($info[0], $info[2]-1, $info[3]+1);
		my @previousInfo = split /\t/, $previousBlock;
		my @postInfo = split /\t/, $postBlock;
		if($info[1] eq 'GAP'){
			($c, $s, $e) = $previousInfo[3] < $postInfo[2] ? ($info[0], $previousInfo[3], $postInfo[2]) : ($info[0], $postInfo[2], $previousInfo[3]);
			$cAlt = $nameTable{$info[0]};
		}
		my $gapString;
		if($c ne $previousChr){
			if($gap{$c}){
				@starts = sort {$a<=>$b} keys %{$gap{$c}};
			}
			elsif($gap{$cAlt}){
				@starts = sort {$a<=>$b} keys %{$gap{$cAlt}};
			}
			$previousChr = $c;
		}
		foreach my $gapS(@starts){
			last if($gapS > $e);
			my $gapE;
			if($gap{$c}){
				$gapE = $gap{$c}{$gapS};
			}
			elsif($gap{$cAlt}){
				$gapE = $gap{$cAlt}{$gapS};
			}
			if(isOverlap($s, $e, $gapS, $gapE)){
				$gapString .= "QN:$gapS-$gapE;";
			}
		}
		print OUT $gapString ? "$l"."$gapString\n" : "$l".".;\n";
	}
}

sub populateNames{
	for my $nameFile($_refNames, $_queryNames){
		open IN, $nameFile;
		while(<IN>){
			chomp;
			my ($name1, $name2) = split /\t/;
			$nameTable{$name1} = $name2;
		}
		close IN;
	}
}