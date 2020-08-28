use strict;
use Getopt::Std;

sub help_mess{
	print "\n", "Analyze Mummer outputs", "\n";
	print "\n";
	print "\t", "-T <Int>", "\t", "Task: ", "\n";
	#print "\t", "", "\t", "1: Sort SV by size", "\n";
	#print "\t", "", "\t", "2: Count Large GAPs", "\n";
	#print "\t", "", "\t", "3: Get Tally PAVs (need -l and -t)", "\n";
	#print "\t", "", "\t", "4: Get PAV sequences from diff file", "\n";
	#print "\t", "", "\t", "5: Count SEQ events", "\n";
	print "\t", "", "\t", "RenameChr: Rename Chromosome. Need -m indicating ID mapping", "\n";
	#print "\t", "", "\t", "8: Get large deletion (compared to reference) from 1 pair-wise alignment", "\n";
	print "\t", "", "\t", "AnnotateGap: Annotate diff and coord file with Gap information. Need -g and -G", "\n";
	print "\t", "", "\t", "NonGapAlignments: Get non-gap events. Need -i and -j to be gap annotated", "\n";
	print "\t", "", "\t", "CountEvents: Count SV events. -i and -j need Gap annotated. Usually need Non-gap alignments", "\n";
	print "\t", "", "\t", "AnnotateGenes: Annotate events file with gene information. Need -r and -q.", "\n";
	print "\t", "", "\t", "AnnotateGapGenes: Annotate genes in filled gaps.", "\n";
	print "\t", "", "\t", "TallyEvents: Identify which event is occuring in which sample. Need file list (-l) of CountEvent outputs.Need -e", "\n";
	print "\t", "", "\t", "GetEventSize: Get event sizes from event files", "\n";
	print "\t", "", "\t", "GetEventSeq: Get event sequences from event files. Need -e and -A", "\n";
	print "\t", "", "\t", "GetSpecificEvents: Output only the specific events for each sample. The input is the tally file from TallyEvents", "\n";
	print "\t", "", "\t", "GetSpecificSequences: Output the sequence associated with specific events. The input is the output files from GetSpecificEvents", "\n";
	print "\t", "", "\t", "AnnotatePAVBlast: Annotate events with blast results. -i is specific events, -l is list of blast outputs.", "\n";
	print "\t", "", "\t", "CheckPAV: Check if PAV based using blast. Need to run GetEventSeq first", "\n";
	print "\t", "", "\t", "CountPAVGenes: Count Genes in PAV regions from event files. If input from specificEvent, need -J", "\n";
	print "\t", "", "\t", "AnnotateDomestication: Annotate whether PAV event falls in domestication regions. need -D", "\n";
	print "\t", "", "\t", "CountSelectedGenes: ", "\n";
	print "\t", "", "\t", "Consolidate: ", "\n";
	print "\t", "", "\t", "TallyConsolidate: ", "\n";
	print "\t", "-i <File>", "\t", "Mummer input file", "\n";
	print "\t", "-l <File>", "\t", "File list", "\n";
	print "\t", "-e [String]", "\t", "Feature to tally", "\n";
	print "\t", "-d [String]", "\t", "Details to record in tally", "\n";
	print "\t", "-t <String>", "\t", "coordinate sorted type (Ref or Query)", "\n";
	print "\t", "-A [String]", "\t", "PAV analysis is for AV or PV", "\n";
	print "\t", "-s [Bool]", "\t", "Get PAV sequence after tallying. Need -f", "\n";
	print "\t", "-S [Bool]", "\t", "Input from specific events", "\n";
	print "\t", "-j [File]", "\t", "Mummer coord file", "\n";
	print "\t", "-J [File]", "\t", "Original Event File", "\n";
	print "\t", "-f [File]", "\t", "Fasta file", "\n";
	print "\t", "-r [File]", "\t", "Ref Gene annotation GFF file", "\n";
	print "\t", "-q [File]", "\t", "Query Gene annotation GFF file", "\n";
	print "\t", "-m [File]", "\t", "ID mapping file", "\n";
	print "\t", "-g [File]", "\t", "Reference Gap info file", "\n";
	print "\t", "-G [File]", "\t", "Query Gap info file", "\n";
	print "\t", "-b [File]", "\t", "Blast output", "\n";
	print "\t", "-D [File]", "\t", "Domestication region file", "\n";
	print "\t", "-p [File]", "\t", "print events after counting", "\n";
	print "\t", "-c [File]", "\t", "chrom name correspondence for Ref", "\n";
	print "\t", "-C	[File]", "\t", "chrom name correspondence for Query", "\n";
	print "\t", "-o [File]", "\t", "output prefix", "\n";
    print "\t", "-h [Bool]", "\t", "Help. Show this and exits.", "\n";
}

## ARGS #############
my %options;
my $arguments = 'i:o:T:j:f:l:t:m:g:G:e:r:q:b:J:D:A:d:c:C:8hspS';
getopts($arguments, \%options);
my $_input = $options{'i'};
my $_output = $options{'o'};
my $_tasks = $options{'T'};
my $_debug = $options{'8'};
my $_help = $options{'h'};
my $_coord_file = $options{'j'};
my $_fasta_file = $options{'f'};
my $_list = $options{'l'};
my $_type = $options{'t'};
my $_get_seq = $options{'s'};
my $_mappingFile = $options{'m'};
my $_rGapFile = $options{'g'};
my $_qGapFile = $options{'G'};
my $_printEvents = $options{'p'};
my $_gffFile = $options{'r'};
my $_qGffFile = $options{'q'};
my $_feature = $options{'e'};
my $_blastFile = $options{'b'};
my $_specificEventInput = $options{'S'};
my $_eventFile = $options{'J'};
my $_domesticationFile = $options{'D'};
my $_PAVType = $options{'A'};
my $_details = $options{'d'};
my $_refNames = $options{'c'};
my $_queryNames = $options{'C'};
%options = ();
if($_help){
	help_mess();
	exit;
}
if(!$_output){
	open OUT, ">&STDOUT";
}
else{
	open OUT, "+>$_output";
}	
if($_get_seq && !$_fasta_file){
	print STDERR "Error: Get sequence requested but no fasta file was provided\n";
	help_mess();
	exit;
}
## Globals #############
my %coord_hash;
my %SV_hash;
my %PAV;
my @keys;
my @files;
my %nameTable;
## Main #############
populateNames();
if($_tasks == 1){
	populate_coords();
	populate_sv();
	#system("grep BRK $_input | sort -nr -k 5 > $_output.BRK.sorted.out");
	system("grep GAP $_input | sort -nr -k 5 > $_output.GAP.sorted.out");
	#sort_GAP();
	if(!$_coord_file){
		print STDERR "Error: Required coord file missing.\n";
		help_mess();
		exit;
	}
	sort_SEQ_INV("SEQ");
	sort_SEQ_INV("INV");
}
elsif($_tasks == 2){
	system("grep GAP $_input | sort -nr -k 5 > $_output.GAP.sorted.out");
	count_GAP();
}
elsif($_tasks == 3){
	tally_pav();
}
elsif($_tasks == 4){
	get_PAV_seqs();
}
elsif($_tasks == 5){
	count_SEQ();
}
elsif($_tasks eq 'CountEvents'){
	CountEvents();
}
elsif($_tasks eq 'RenameChr'){
	RenameChrV2();
}
elsif($_tasks eq 'AnnotateGap'){
	if(!$_coord_file){
		print STDERR "Error: Required coord file missing.\n";
		help_mess();
		exit;
	}
	if(!$_rGapFile || !$_qGapFile){
		print STDERR "Error: Required Gap file missing.\n";
		help_mess();
		exit;
	}
	AnnotateRGap();
	AnnotateQGap();
	my($prefix, $suffix) = $_coord_file =~ /(.*)\.(\w+)/;
	system("cp $prefix.rGapAnno.$suffix $prefix.gapAnno.$suffix");
	($prefix, $suffix) = $_input =~ /(.*)\.(\w+)/;
	#system("rm $prefix.rGapAnno.$suffix");
}
elsif($_tasks eq 'NonGapAlignments'){
	if(!$_coord_file){
		print STDERR "Error: Required coord file missing.\n";
		help_mess();
		exit;
	}
	NonGapAlignmentsDiff();
	NonGapAlignmentsCoord();
}
elsif($_tasks eq 'AnnotateGenes'){
	if(!$_gffFile){
		print STDERR "Error: Required Gff file missing.\n";
		help_mess();
		exit;
	}
	AnnotateGenes();
}
elsif($_tasks eq 'TallyEvents'){
	die "Error: Required list file missing\n" if(!$_list);
	die "Error: Required feature to tally missing\n" if(!$_feature);
	die "Error: Specificing PV or AV required\n" if(!$_PAVType);
	TallyEvents();
}
elsif($_tasks eq 'Consolidate'){
	if(!$_list || !$_feature){
		print STDERR "Error: Required list file or feature missing.\n";
		help_mess();
		exit;
	}
	consolidate();
}
elsif($_tasks eq 'TallyConsolidate'){
	TallyConsolidate();
}
elsif($_tasks eq 'GetSpecificEvents'){
	GetSpecificEvents();
}
elsif($_tasks eq 'GetSpecificSequences'){
	GetSpecificSequences();
}
elsif($_tasks eq 'AnnotatePAVBlast'){
	AnnotatePAVBlast();
}
elsif($_tasks eq 'CheckPAV'){
	CheckPAV($_PAVType);
}
elsif($_tasks eq 'GetEventSeq'){
	if(!$_feature || !$_PAVType){
		print STDERR "Error: Required feature missing.\n";
		help_mess();
		exit;
	}
	GetEventSeq($_PAVType);
}
elsif($_tasks eq 'GetEventSize'){
	if($_specificEventInput){
		GetEventSizeCoord();
	}
	else{
		GetEventSize();
	}
}
elsif($_tasks eq 'CountPAVGenes'){
	if($_specificEventInput){
		CountPAVGenesCoord();
	}
	else{
		CountPAVGenes();
	}
}
elsif($_tasks eq 'AnnotateDomestication'){
	AnnotateDomestication();
}
elsif($_tasks eq 'CountSelectedGenes'){
	if($_specificEventInput){
		CountSelectedGenesCoord();
	}
	else{
		CountSelectedGenes();
	}
}
elsif($_tasks eq 'AnnotateGapGenes'){
	if(!$_input){
		print STDERR "Error: Required filled gap input file missing.\n";
		help_mess();
		exit;
	}
	if(!$_gffFile){
		print STDERR "Error: Required ref gff file missing.\n";
		help_mess();
		exit;
	}
	AnnotateGapGenes();
}
else{
	print STDERR "Error: Unknown Task.\n";
		help_mess();
		exit;
}
## Subroutines #############
sub CountSelectedGenesCoord{
	open IN, $_input;
	my %hash;
	my %genes;
	my $totalGenes;
	while(my $l = <IN>){
		chomp $l;
		my ($refCoord, $queryCoord) = split /\t/, $l;
		$hash{$queryCoord} = 1;
	}	
	close IN;
	
	open IN, $_eventFile;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;		
		if($info[11] eq 'PAV'){
			if($info[12] ne '-' || $info[13] ne '-'){
				if($info[10] ne '-'){
					if($hash{"$info[4]:$info[5]-$info[6]"}){
						my @genes = split /;/, $info[10];
						foreach my $g(@genes){
							$genes{$g} = 1;
						}
						$totalGenes += $#genes+1;
						print STDERR $l, "\n";
					}					
				}
			}
		}
	}	
	print "Total Selected Genes in Specific PAV = $totalGenes\n";
	print "Genes:\n", join("\n", keys %genes), "\n";
	
}

sub CountSelectedGenes{
	open IN, $_input;
	my $totalGenes;
	my %genes;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		if($info[11] eq 'PAV'){
			if($info[12] ne '-' || $info[13] ne '-'){
				if($info[10] ne '-'){
					my @genes = split /;/, $info[10];
					$totalGenes += $#genes+1;
					foreach my $g(@genes){
						$genes{$g} = "$info[1]:$info[2]-$info[3]";
					}
				}
			}
		}
	}	
	print "Total Selected Genes in PAV = $totalGenes\n";
	foreach my $k(keys %genes){
		print $k, "\t", $genes{$k}, "\n";
	}
}

sub AnnotateDomestication{
	my $type;
	if($_domesticationFile =~ /domestication\.\d+\.region/){
		$type = 'domestication';
	}
	elsif($_domesticationFile =~ /improvement\.\d+\.region/){
		$type = 'improvement';
	}
	open IN, $_domesticationFile;
	my %hash;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my ($refChrom, $refStart, $refEnd) = $info[3] < $info[4] ? ($info[2], $info[3], $info[4]) : ($info[2], $info[4], $info[3]);
		$hash{$refChrom}{$refStart} = $refEnd;
	}
	close IN;
	
	open IN, $_input;
	my $previousChr;
	my @starts;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		next if($info[0] ne 'Indel');
		my $overlapStr = '-';
		my ($refChrom, $refStart, $refEnd) = $info[2] < $info[3] ? ($info[1], $info[2], $info[3]) : ($info[1], $info[3], $info[2]);
		my ($pavC, $pavS, $pavE) = $info[5] < $info[6] ? ($info[4], $info[5], $info[6]) : ($info[4], $info[6], $info[5]);
		#my ($currChrom, $currS, $currE) = ($pavC, $pavS, $pavE);
		my ($currChrom, $currS, $currE) = ($refChrom, $refStart, $refEnd);
		if($currChrom ne $previousChr){
			@starts = keys %{$hash{$currChrom}};
			$previousChr = $currChrom;
		}
		foreach my $s(@starts){
			if(isOverlap($currS, $currE, $s, $hash{$currChrom}{$s})){
				$overlapStr = $type;
				last;
			}
		}
		print $l, "\t", $overlapStr, "\n";
	}
}

sub CountPAVGenesCoord{
	open IN, $_input;
	my %hash;
	my $totalGenes;
	while(my $l = <IN>){
		chomp $l;
		my ($refCoord, $queryCoord) = split /\t/, $l;
		$hash{$queryCoord} = 1;
	}	
	close IN;
	
	open IN, $_eventFile;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		if($info[10] ne '-'){
			if($hash{"$info[4]:$info[5]-$info[6]"}){
				my @genes = split /;/, $info[10];
				$totalGenes += $#genes+1;
			}
		}
	}	
	print "Total Genes in Specific PAV = $totalGenes\n";
}

sub CountPAVGenes{
	open IN, $_input;
	my $totalGenes;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		if($info[10] ne '-'){
			my @genes = split /;/, $info[10];
			$totalGenes += $#genes+1;
		}
	}
	print "Total Genes in PAV = $totalGenes\n";
}

sub GetEventSizeCoord{
	open IN, $_input;
	my ($totalRefSize, $totalQuerySize);
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my ($rc, $rs, $re) = $info[0] =~ /(.*):(\d+)-(\d+)/;
		my ($qc, $qs, $qe) = $info[1] =~ /(.*):(\d+)-(\d+)/;
		my $rSize = abs($re - $rs + 1);
		my $qSize = abs($qe - $qs + 1);
		$totalRefSize += $rSize;
		$totalQuerySize += $qSize;
	}
	print "Total Ref Size = $totalRefSize\nTotal Query Size = $totalQuerySize\n";
}

sub GetEventSize{
	open IN, $_input;
	my ($totalRefSize, $totalQuerySize);
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my $rSize = abs($info[3] - $info[2] + 1);
		my $qSize = abs($info[6] - $info[5] + 1);
		$totalRefSize += $rSize;
		$totalQuerySize += $qSize;
	}
	print "Total Ref Size = $totalRefSize\nTotal Query Size = $totalQuerySize\n";
}

sub GetEventSeq{
	my $type = shift;
	open IN, $_input;
	open OUT, "+>$_input.$type.tmp";
	while(my $l = <IN>){
		my @info = split /\t/, $l;
		if($type eq 'AV'){
			if($info[3]-$info[2]+1 < 100 ||
				$info[0] ne $_feature){
				next;
			}
			print OUT "$info[1]:$info[2]-$info[3]", "\n";
		}
		elsif($type eq 'PV'){
			if($info[6]-$info[5]+1 < 100 ||
				$info[0] ne $_feature){
				next;
			}
			print OUT "$info[4]:$info[5]-$info[6]", "\n";
		}
	}
	#system("perl ~/scripts/fasta_tools.pl -T 4 -i $_fasta_file -j $_input.$type.tmp -o $_input.$type.tmp.fasta");
}

sub CheckPAV{
	my $type = shift;
	my %hash;
	open IN, $_blastFile;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my $coverage = ($info[7] - $info[6] + 1) / $info[12] * 100;
		if($info[2] > 90 && $coverage > 50){
			$hash{$info[0]} = '-';
		}
		else{
			$hash{$info[0]} = 'PAV';
		}
	}	
	close IN;
	
	open IN, $_input;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my $coord;
		if($type eq 'AV'){
			if($info[3]-$info[2]+1 < 100 || $info[0] ne 'Indel'){
				print "$l\t-\n";
				next;
			}
			$coord = "$info[1]:$info[2]-$info[3]";
		}
		elsif($type eq 'PV'){
			if($info[6]-$info[5]+1 < 100 || $info[0] ne 'Indel'){
				print "$l\t-\n";
				next;
			}
			$coord = "$info[4]:$info[5]-$info[6]";
		}
		
		if($hash{$coord}){
			print $l, "\t", $hash{$coord}, "\n";
		}
		else{
			print $l, "\t-\n";
		}
	}
	close IN;
}

sub AnnotatePAVBlast{
	my %hash;
	open IN, $_input;
	while(my $l = <IN>){
		chomp $l;
		my($ref, $query) = split /\t/, $l;
		$hash{$query} = 'PAV';
	}
	close IN;
	open LI, $_list;
	while(my $file = <LI>){
		chomp $file;
		open IN, $file;
		while(my $l = <IN>){
			chomp $l;
			my @info = split /\t/, $l;
			
		}
	}
}
sub GetSpecificEvents{
	open IN, $_input;
	chomp(my $header = <IN>);
	my @header = split /\t/, $header;
	my @events;

	while(my $l = <IN>){
		chomp $l;
		my ($count, $s);
		my @info = split /\t/, $l;
		for(my $i = 1; $i < @header; $i++){
			if($info[$i] ne '.'){
				$count++;
				$s = $i;
			}
		}
		if($count == 1){
			$events[$s] .= "$info[0]\t$info[$s]\n";
		}	
	}
	for(my $i = 1; $i < @header; $i++){
		open OUT, "+>$_input.$i.specificEvent";
		print OUT $events[$i], "\n";
	}
}

sub consolidate{
	open LI, $_list;
	my $index = 1;
	while(my $file = <LI>){
		chomp $file;
		push @files, $file;
		print STDERR "Tallying $file\n";
		open IN, $file;
		my $previousChr;
		while(my $l = <IN>){
			chomp $l;
			my @info = split /\t/, $l;
			next if($info[0] ne $_feature);
			next if($info[7] < 100 && $info[8] < 100);
			my ($c, $s, $e) = $info[2] < $info[3] ? ($info[1], $info[2], $info[3]) : ($info[1], $info[3], $info[2]);
			my $hasOverlap = 0;
			my $details;
			if($_details eq 'size'){
				$details = max($info[7], $info[8]);
			}
			elsif($_details eq 'coord'){
				$details = "$info[4]:$info[5]-$info[6]";
			}
			elsif($_details eq 'genes'){
				$details  = "$info[9]|$info[10]";
			}
			else{
				$details = "$info[4]:$info[5]-$info[6]";
			}
			if($index == 1){
				$PAV{$c}{$s}[0] = $e;
				$PAV{$c}{$s}[1] = $details;
			}
			else{
				if($PAV{$c}{$s}){
					$PAV{$c}{$s}[0] = $e if($e > $PAV{$c}{$s}[0]);
					if($_details eq 'size'){
						$PAV{$c}{$s}[1] = $details if $details > $PAV{$c}{$s}[1];
					}
					else{
						$PAV{$c}{$s}[1] .= $details;
					}
					
				}
				else{
					if($c ne $previousChr){
						@keys = sort {$a<=>$b} keys %{$PAV{$c}};
						$previousChr = $c;
					}				
					foreach my $eventS(@keys){
						last if($eventS > $e);
						if(isOverlap($s, $e, $eventS, $PAV{$c}{$eventS}[0])){
							$PAV{$c}{$eventS}[0] = $e if($e > $PAV{$c}{$eventS}[0]);
							if($_details eq 'size'){
								$PAV{$c}{$s}[1] = $details if $details > $PAV{$c}{$s}[1];
							}
							else{
								$PAV{$c}{$eventS}[1] .= $details;
							}
							$hasOverlap = 1;
						}
					}
					if(!$hasOverlap){
						$PAV{$c}{$s}[0] = $e;
						$PAV{$c}{$s}[1] = $details;
					}	
				}
			}
		}
		$index++;
	}

	#Print tally table
	my @specific;
	#print "Reference Coord\t", join("\t", @files), "\n";
	foreach my $c(sort keys %PAV){
		foreach my $s(sort {$a <=> $b} keys %{$PAV{$c}}){
			my ($count, $sample, $genes);
			print "$_feature\t$c:$s-$PAV{$c}{$s}[0]\t$PAV{$c}{$s}[1]\n";
		}
	}
}

sub TallyConsolidate{
	open LI, $_list;
	my $index = 1;
	while(my $file = <LI>){
		chomp $file;
		push @files, $file;
		print STDERR "Tallying $file\n";
		open IN, $file;
		my $previousChr;
		while(my $l = <IN>){
			chomp $l;
			my @info = split /\t/, $l;
			my ($c, $s, $e) = $info[0] =~ /(\S+):(\d+)-(\d+)/;
			my $hasOverlap = 0;
			my $details = $info[1];
			if($index == 1){
				$PAV{$c}{$s}[0] = $e;
				$PAV{$c}{$s}[$index] = $details;
			}
			else{
				if($PAV{$c}{$s}){
					$PAV{$c}{$s}[0] = $e if($e > $PAV{$c}{$s}[0]);
					$PAV{$c}{$s}[$index] = $details if $details > $PAV{$c}{$s}[$index];
				}
				else{
					if($c ne $previousChr){
						@keys = sort {$a<=>$b} keys %{$PAV{$c}};
						$previousChr = $c;
					}				
					foreach my $eventS(@keys){
						last if($eventS > $e);
						if(isOverlap($s, $e, $eventS, $PAV{$c}{$eventS}[0])){
							$PAV{$c}{$eventS}[0] = $e if($e > $PAV{$c}{$eventS}[0]);
							$PAV{$c}{$eventS}[$index] = $details if $details > $PAV{$c}{$eventS}[$index];
							$hasOverlap = 1;
						}
					}
					if(!$hasOverlap){
						$PAV{$c}{$s}[0] = $e;
						$PAV{$c}{$s}[$index] = $details;
					}	
				}
			}
			#if(!$PAV{$c}{$s}[0]){
			#	print STDERR "End missing for $c:$s\n";
			#	<>;
			#}
		}
		$index++;
	}
	#Print tally table
	my @specific;
	print "Reference Coord\t", join("\t", @files), "\n";
	foreach my $c(sort keys %PAV){
		foreach my $s(sort {$a <=> $b} keys %{$PAV{$c}}){
			my ($count, $sample, $genes);
			print "$c:$s-$PAV{$c}{$s}[0]\t";
			for(my $i = 1; $i <= @files; $i++){
				if($PAV{$c}{$s}[$i]){
					$count++;
					$genes = $PAV{$c}{$s}[$i] if($_details eq 'genes');
					print $PAV{$c}{$s}[$i], "\t";
					$sample = $i;
				}
				else{
					print "0\t";
				}
			}
			$genes = '-' unless $genes;
			print $count, "\t", $genes, "\n";
			if($count == 1){
				$specific[$sample-1]++;
			}			
		}
	}
}


sub TallyEvents{
	open LI, $_list;
	my $index = 1;
	while(my $file = <LI>){
		chomp $file;
		push @files, $file;
		print STDERR "Tallying $file\n";
		open IN, $file;
		my $previousChr;
		while(my $l = <IN>){
			chomp $l;
			#if($file =~ /08/){
			#	print STDERR "line = $l\n";
			#}
			my @info = split /\t/, $l;
			next if($info[0] ne $_feature);
			next if($info[7] < 100 && $info[8] < 100);
			#next if($info[1] ne $info[4]);
			my ($c, $s, $e) = $info[2] < $info[3] ? ($info[1], $info[2], $info[3]) : ($info[1], $info[3], $info[2]);
			my $hasOverlap = 0;
			my $details;
			#my $details = "\[$info[9]|$info[10]\]";
			#my $details = "$info[4]:$info[5]-$info[6]";
			#my $details  = "$info[10]";
			if($_details eq 'size'){
				#$details = max($info[7], $info[8]);
				if($_PAVType eq 'AV'){
					$details = $info[7];
				}
				elsif($_PAVType eq 'PV'){
					$details = $info[8];
				}
			}
			elsif($_details eq 'coord'){
				$details = "$info[4]:$info[5]-$info[6]";
			}
			elsif($_details eq 'genes'){
				$details  = "$info[9]|$info[10]";
			}
			else{
				$details = "$info[4]:$info[5]-$info[6]";
			}
			if($index == 1){
				$PAV{$c}{$s}[0] = $e;
				$PAV{$c}{$s}[$index] = $details;
			}
			else{
				if($PAV{$c}{$s}){
					$PAV{$c}{$s}[0] = $e if($e > $PAV{$c}{$s}[0]);
					if($_details eq 'size'){
						$PAV{$c}{$s}[$index] = $details if $details > $PAV{$c}{$s}[$index];
					}
					else{
						$PAV{$c}{$s}[$index] .= "$details";
					}
					
				}
				else{
					if($c ne $previousChr){
						@keys = sort {$a<=>$b} keys %{$PAV{$c}};
						$previousChr = $c;
					}				
					foreach my $eventS(@keys){
						last if($eventS > $e);
						if(isOverlap($s, $e, $eventS, $PAV{$c}{$eventS}[0])){
							$PAV{$c}{$eventS}[0] = $e if($e > $PAV{$c}{$eventS}[0]);
							if($_details eq 'size'){
								$PAV{$c}{$eventS}[$index] = $details if $details > $PAV{$c}{$eventS}[$index];
							}
							else{
								$PAV{$c}{$eventS}[$index] .= "$details";
							}
							$hasOverlap = 1;
						}
					}
					if(!$hasOverlap){
						$PAV{$c}{$s}[0] = $e;
						$PAV{$c}{$s}[$index] = $details;
					}	
				}
			}
		}
		$index++;
	}
	
	while(check_overlap()){
		print STDERR "Cleaning up overlaps in hash\n";
	}
	print STDERR "No more overlaps found\n";

	#Print tally table
	my @specific;
	print "Reference Coord\t", join("\t", @files), "\n";
	foreach my $c(sort keys %PAV){
		foreach my $s(sort {$a <=> $b} keys %{$PAV{$c}}){
			my ($count, $sample, $genes);
			print "$c:$s-$PAV{$c}{$s}[0]\t";
			for(my $i = 1; $i <= @files; $i++){
				if($PAV{$c}{$s}[$i]){
					$count++;
					$genes = $PAV{$c}{$s}[$i] if($_details eq 'genes');
					print $PAV{$c}{$s}[$i], "\t";
					$sample = $i;
				}
				else{
					if($_details eq 'size'){
						print "0\t";
					}
					else{
						print ".\t";
					}
				}
			}
			$genes = '-' unless $genes;
			#print $count, "\t", $genes, "\n";
			print "\n";
			if($count == 1){
				$specific[$sample-1]++;
			}			
		}
	}
	for(my $i = 0; $i < @files; $i++){
		#print STDERR $files[$i], "\t", $specific[$i], "\n";
	}
}

sub AnnotateGenes{
	my %hash;
	open IN, $_gffFile;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		$info[0] =~ s/_rc//;
		next unless($info[2] eq 'CDS');
		my ($name) = $info[8] =~ /Parent=(.*?);/;
		#print STDERR "$info[0], $info[3], $info[4], $name\n";
		$hash{$info[0]}{$info[3]}[0] = $info[4];
		$hash{$info[0]}{$info[3]}[1] = $name;
	}
	close IN;
	
	my %queryHash;
	open IN, $_qGffFile;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		$info[0] =~ s/_rc//;
		next unless($info[2] eq 'CDS');
		my ($name) = $info[8] =~ /Parent=(.*?);/;
		#print STDERR "$info[0], $info[3], $info[4], $name\n";
		$queryHash{$info[0]}{$info[3]}[0] = $info[4];
		$queryHash{$info[0]}{$info[3]}[1] = $name;
	}
	close IN;	
	
	open IN, $_input;
	my ($previousChr, $previousQChr);
	my (@starts, @qStarts);
	my $gapRGenes;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my (%rGenes, %qGenes);
		my ($c, $s, $e) = $info[2] < $info[3] ? ($info[1], $info[2], $info[3]) : ($info[1], $info[3], $info[2]);
		my ($qC, $qS, $qE) = $info[5] < $info[6] ? ($info[4],$info[5],$info[6]) : ($info[4], $info[6], $info[5]);
		if($previousChr ne $c){
			@starts = sort {$a<=>$b} keys %{$hash{$c}};
			$previousChr = $c;
			#print STDERR "Getting new array for $c, array length = $#starts\n";
		}
		if($previousQChr ne $qC){
			@qStarts = sort {$a<=>$b} keys %{$queryHash{$qC}};
			$previousQChr = $qC;
		}
		foreach my $refS(@starts){
			#print STDERR "Compare $s-$e and $refS-$hash{$c}{$refS}[0]\n";
			last if($refS > $e);
			if(isOverlap($s, $e, $refS, $hash{$c}{$refS}[0])){
				$rGenes{$hash{$c}{$refS}[1]} = 1;
			}
		}
		foreach my $queryS(@qStarts){
			last if($queryS > $qE);
			if(isOverlap($qS, $qE, $queryS, $queryHash{$qC}{$queryS}[0])){
				$qGenes{$queryHash{$qC}{$queryS}[1]} = 1;
			}
		}
		my @rGenes = keys %rGenes;
		my @qGenes = keys %qGenes;
		if(@rGenes){
			print $l, "\t", join(";", @rGenes), "\t";
			if($info[0] eq 'GapR'){
				#print STDERR "Length 1 = ", $#genes + 1, "\n";
				#print STDERR "Length 2 = ", length @genes, "\n";
				$gapRGenes += length @rGenes;
			}
		}
		else{
			print $l, "\t-\t"; 
		}
		if(@qGenes){
			print join(";", @qGenes), "\n";
		}
		else{
			print "-\n";
		}
	}
	print STDERR "$gapRGenes Genes in GapR regions\n";
}


sub AnnotateGapGenes{
	my %hash;
	open IN, $_gffFile;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		$info[0] =~ s/_rc//;
		next unless($info[2] eq 'CDS');
		my ($name) = $info[8] =~ /Parent=(.*?);/;
		#print STDERR "$info[0], $info[3], $info[4], $name\n";
		$hash{$info[0]}{$info[3]}[0] = $info[4];
		$hash{$info[0]}{$info[3]}[1] = $name;
	}
	close IN;
	
	open IN, $_input;
	my ($previousChr, $previousQChr);
	my (@starts, @qStarts);
	my $gapRGenes;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my (%rGenes);
		my ($c, $s, $e) = ($info[0], $info[2], $info[3]);
		
		if($previousChr ne $c){
			@starts = sort {$a<=>$b} keys %{$hash{$c}};
			$previousChr = $c;
			#print STDERR "Getting new array for $c, array length = $#starts\n";
		}
		foreach my $refS(@starts){
			last if($refS > $e);
			next if($hash{$c}{$refS}[0] < $s);
			if(isOverlap($s, $e, $refS, $hash{$c}{$refS}[0])){
				$rGenes{$hash{$c}{$refS}[1]} = 1;
			}
		}
		my @rGenes = keys %rGenes;
		if(@rGenes){
			print $l, "\t", join(";", @rGenes), "\n";
		}
		else{
			print $l, "\t-\n"; 
		}
	}
	
}


sub NonGapAlignmentsDiff{
	open IN, $_input;
	my($prefix, $suffix) = $_input =~ /(.*)\.(\w+)/;
	open OUT, "+>$prefix.nonGap.$suffix";
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;		
		my $previous = $info[2] - 1;
		my $post = $info[3] + 1;
		chomp(my $previous_block = `grep -P '\t$previous\t.*$info[0]\t' $_coord_file`);
		chomp(my $post_block = `grep -P '^$post\t.*$info[0]\t' $_coord_file`);
		my @previous_info = split /\t/, $previous_block;
		my @post_info = split /\t/, $post_block;
		if($info[1] eq 'GAP'){
			print OUT "$l\n"; #prints all gaps for counting stats later in count Events.
		}
		elsif($info[1] eq 'SEQ'){
			if($info[5] eq $info[0] && $info[6] ne $info[0]){
				#Opening SEQ
				print OUT "$l\n" unless $post_block =~ /RN/;
			}
			elsif($info[5] ne $info[0] && $info[6] eq $info[0]){
				#closing SEQ
				print OUT "$l\n" unless $previous_block =~ /RN/;
			}
			else{
				print OUT "$l\n" unless($previous_block =~ /RN/ || $post_block =~ /RN/);
			}
		}
		elsif($info[1] eq 'JMP'){
			print OUT "$l\n" unless($previous_block =~ /RN/ || $post_block =~ /RN/);
		}
		elsif($info[1] eq 'INV'){
			if($previous_info[11] eq $previous_info[12] && $post_info[11] ne $post_info[12]){
				#Opening INV
				print OUT "$l\n" unless $post_block =~ /RN/;
			}
			elsif($previous_info[11] ne $previous_info[12] && $post_info[11] eq $post_info[12]){
				#Closing INV
				print OUT "$l\n" unless $previous_block =~ /RN/;
			}
			else{
				print OUT "$l\n" unless($previous_block =~ /RN/ || $post_block =~ /RN/);
			}
		}
		elsif($info[1] eq 'BRK'){
			#print "$l\n" unless($l =~ /RN/);
			print OUT "$l\n";
		}
	}
	close IN;
	close OUT;
}

sub NonGapAlignmentsCoord{
	open IN, $_coord_file;
	my($prefix, $suffix) = $_coord_file =~ /(.*)\.(\w+)/;
	open OUT, "+>$prefix.nonGap.$suffix";
	while(my $l = <IN>){
		chomp $l;
		next if($l =~ /RN/);
		print OUT $l, "\n";
	}		
	close IN;
	close OUT;
}

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
	
	#Annotate diff file (Only sort by ref coordinate)
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
	open IN, $_qGapFile or die "Can't find qGapFile $_qGapFile";
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
=cut		
		elsif($info[1] eq 'INV'){
			if($previousInfo[13] ne $previousInfo[14] || $postInfo[13] ne $postInfo[14]){
				next;
			}
			elsif($previousInfo[11] eq $previousInfo[12]){
				#Opening INV
				($invC, $invE) = ($previousInfo[14], $previousInfo[3]);
			}	
			elsif($postInfo[11] eq $postInfo[12]){
				#Closing INV
				($c, $s, $e) = ($invC, $postInfo[2], $invE);
			}
		}
		next unless ($c && $s && $e);
=cut
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

sub getBlocks{
	my ($chr, $previous, $post) = @_;
	chomp(my $prB = `grep -P '\t$previous\t.*$chr\t' $_coord_file`);
	chomp(my $poB = `grep -P '^$post\t.*$chr\t' $_coord_file`);
	print STDERR "WARNING: Cannot find previous Block using '\t$previous\t.*$chr\t'\n" if !$prB;
	print STDERR "WARNING: Cannot find post block using '^$post\t.*$chr\t'\n" if !$poB;
	return ($prB, $poB);
}

sub isOverlap{
	my ($s1, $e1, $s2, $e2) = @_;
	if(max($s1, $s2) <= min($e1, $e2)){
		return 1;
	}
	return 0;
}

sub RenameChr{
	my %name;
	open IN, $_mappingFile;
	while(<IN>){
		chomp;
		my ($id1, $id2) = split /\s+/;
		$name{$id1} = $id2;
	}
	close IN;
	my @names = keys %name;
	open LI, $_list;
	while(my $file = <LI>){
		chomp $file;
		my ($prefix, $suffix) = $file =~ /(.*)\.(\w+)/;
		open IN, $file;
		open OUT, "+>$prefix.rename.$suffix";
		while(my $l = <IN>){
			chomp $l;
			foreach my $n(@names){
				#print STDERR "n1 = $n, n2 = $name{$n}\n";
				$l =~ s/$n/$name{$n}/;
			}
			print OUT $l, "\n";
		}
		close IN;
	}
	close LI;
}

sub RenameChrV2{
	my %name;
	open IN, $_mappingFile;
	while(<IN>){
		chomp;
		my ($id1, $id2) = split /\s+/;
		$name{$id1} = $id2;
	}
	close IN;
	my @names = keys %name;
	open LI, $_list;
	while(my $file = <LI>){
		chomp $file;
		my ($prefix, $suffix) = $file =~ /(.*)\.(\w+)/;
		open IN, $file;
		open OUT, "+>$prefix.rename.$suffix";
		read IN, my $content, -s IN;
		foreach my $n(@names){
			print STDERR "Substituting $n for $name{$n}\n";
			$content =~ s/$n/$name{$n}/gm;
		}
		print OUT $content, "\n";
		close IN;
	}
	close LI;
}

sub get_indel{
	open IN, $_input;
	open PV, "+>$_input.PV";
	open AV, "+>$_input.AV";
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		next if($info[1] ne 'GAP');
		my ($s, $e) = $info[2] < $info[3] ? ($info[2], $info[3]) : ($info[3], $info[2]);
		if($info[4] >= 100){	
			print AV "$info[0]\t$s\t$e\t", $e-$s+1, "\n";
		}
		if($info[4] < 0 && $info[5] > 100){
			print PV "$info[0]\t$s\t$e\n";
		}
	}
}

sub CountEvents{
	my ($invtrans, $previousChr);
	my ($totalIndel, $repeatIndel, $gap, $totalLength, $validLength);
	my ($gapBoth, $gapR, $gapQ);
	my ($translocation, $translocationGap, $translocationLength);
	my ($seqStart, $seqEnd, $seqGap);
	my ($inversion, $inversionGap, $inversionLength);
	my ($invStart, $invEnd, $invGap);
	open IN, $_input;
	while(my $l = <IN>){
		chomp $l;
		#next if($l =~ /[sS]caffold/ || $l =~ /tig/);
		my @info = split /\t/, $l;
		my $type;
		if($info[0] ne $previousChr){
			($seqStart, $seqEnd, $seqGap) = ();
			($invStart, $invEnd, $invGap) = ();
			$previousChr = $info[0];
		}
		if($info[1] eq 'GAP'){
			my ($previousBlock, $postBlock) = getBlocks($info[0], $info[2]-1, $info[3]+1);
			next unless($previousBlock && $postBlock);
			my @previousInfo = split /\t/, $previousBlock;
			my @postInfo = split /\t/, $postBlock;
			next if($nameTable{$info[0]} ne $nameTable{$previousInfo[14]} || $nameTable{$info[0]} ne $nameTable{$postInfo[14]});
			#next if($info[0] ne $previousInfo[14] || $info[0] ne $postInfo[14]); # next if it is translocation. Need to fix if naming covention do not match
			$previousInfo[3]++;
			$postInfo[2]--;
			$totalIndel++;
			$totalLength += abs($info[6]);
			if($l =~ /RN.*QN/){
				$gap++;
				$gapBoth++;
				$type = 'GapBoth';
			}
			elsif($l =~ /RN/){
				$gap++;
				$gapR++;
				$type = 'GapR';
			}
			elsif($l =~ /QN/){
				$gap++;
				$gapQ++;
				$type = 'GapQ';
			}
			elsif($info[4] < 0 && $info[5] < 0){
				$repeatIndel++;
				$type = 'RepetitiveIndel';
			}
			else{
				$type = 'Indel';
				$validLength += abs($info[6]);
				if($info[4] > 50 || $info[5] > 50){
					print "$type\t$info[0]\t$info[2]\t$info[3]\t$previousInfo[14]\t$previousInfo[3]\t$postInfo[2]\t$info[4]\t$info[5]\n";
				}
			}
		}	
		elsif($info[1] eq 'SEQ'){
			if($info[5] eq $info[0] && $info[6] ne $info[0]){
				#Opening SEQ
				$seqStart = $info[2] < $info[3] ? $info[2] : $info[3];
				$seqGap = 1 if $l =~ /RN/;
			}
			elsif($info[5] ne $info[0] && $info[6] eq $info[0]){
				#closing SEQ
				if($seqStart){
					$seqEnd = $info[2] < $info[3] ? $info[3] : $info[2];
					$seqGap = 1 if $l =~ /RN/;
					$translocationLength += $seqEnd - $seqStart + 1;
					$seqGap ? $translocationGap++ : $translocation++;
					#print "Translocation\t$info[0]\t$seqStart\t$seqEnd\n";
					($seqStart, $seqEnd, $seqGap) = ();
				}
			}	
		}
=cut
		elsif($info[1] eq 'INV'){
			my $previous = $info[2] - 1;
			my $post = $info[3] + 1;
			chomp(my $previous_block = `grep -P '\t$previous\t.*$info[0]\t' $_coord_file`);
			chomp(my $post_block = `grep -P '^$post\t.*$info[0]\t' $_coord_file`);
			my @previous_info = split /\t/, $previous_block;
			my @post_info = split /\t/, $post_block;
			if($previous_info[11] eq $previous_info[12] && $post_info[11] ne $post_info[12]){
				#Opening INV
				$invStart = $info[2] < $info[3] ? $info[2] : $info[3];
				$invGap = 1 if $l =~ /RN/;
				$invEnd = ();
			}
			elsif($previous_info[11] ne $previous_info[12] && $post_info[11] eq $post_info[12]){
				#Closing INV
				if($invStart){
					$invGap = 1 if $l =~ /RN/;
					$invEnd = $info[2] < $info[3] ? $info[3] : $info[2];
					$inversionLength += $invEnd - $invStart + 1;
					$invGap ? $inversionGap++ : $inversion++;
					print "Inversion\t$info[0]\t$invStart\t$invEnd\n";
					($invStart, $invEnd, $invGap) = ();
				}
			}
		}
=cut
	}	
	close IN;
	#Print inversion events
	($inversion, $inversionGap, $inversionLength) = printInversion();
	#Print translocation events
	($translocation, $translocationGap, $translocationLength) = printTranslocation();

	print STDERR "$_input and $_coord_file has:\n";
	print STDERR "$totalIndel Total Indels\n$repeatIndel Indels in Repeat region\n$gap Gapped Indels\nTotal Indel length: $totalLength\n";
	print STDERR $totalIndel - $repeatIndel - $gap, " Valid Indels\nValid Indel length: $validLength\n";
	print STDERR "Of the $gap Gap events:\n$gapBoth have gap in both genomes\n$gapR have gap in Ref genome\n$gapQ have gap in Query genome\n";
	print STDERR "----------\n";
	print STDERR $translocation + $translocationGap, " Total translocations\n$translocation Valid Translocations\n$translocationGap Translocations with Gapped boundaries\nTotal translocation length: $translocationLength\n";
	print STDERR "----------\n";
	print STDERR $inversion + $inversionGap, " Total inversions\n$inversion Inversions\n$inversionGap Inversions with Gapped Boundaries\nTotal inversion length: $inversionLength\n";
	#print "$invtrans Inversion&Translocations\n";
	print STDERR "----------\n";
	print STDERR "Total events = ", $totalIndel + $translocation + $translocationGap + $inversion + $inversionGap, "\n";
	print STDERR "Total valid events = ", $totalIndel - $repeatIndel - $gap + $translocation + $inversion, "\n";
	print STDERR "Total gapped events = ", $gap + $translocationGap + $inversionGap, "\n";
	print STDERR "Total event length = ", $totalLength + $translocationLength + $inversionLength, "\n";
	print STDERR "Total valid event length = ", $validLength + $translocationLength + $inversionLength, "\n";
}

sub getDiffs{
	my ($chr, $previous, $post) = @_;
	chomp(my $prB = `grep -P '^$chr\t.*\t$previous\t' $_input`);
	chomp(my $poB = `grep -P '^$chr\t.*\t$post\t' $_input`);
	return ($prB, $poB);
}


sub printInversion{
	open IN, $_coord_file;
	my $previous = 0;
	my ($inv, $invS, $invE, $invLength, $previousChr, $invGap, $type);
	my ($altC, $altS, $altE);
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		next if($nameTable{$info[13]} ne $nameTable{$info[14]});
		#next if($info[13] ne $info[14]);
		#next if($l =~ /[sS]caffold/ || $l =~ /tig/);
		if($info[11] ne $info[12]){
			if($info[13] ne $previousChr){
				#Different Chromosome
				if($previous){
					#End current INV event
					$invLength += $invE - $invS + 1;
					my ($previousDiff, $postDiff) = getDiffs($previousChr, $invS-1, $invE+1);
					if($previousDiff =~ /[RQ]N/ || $postDiff =~ /[RQ]N/){
						$invGap++;
						$type = 'InversionGap';
					}
					else{
						$inv++;
						$type = 'Inversion';
						print "$type\t$previousChr\t$invS\t$invE\t$altC\t$altS\t$altE\t", $invE-$invS+1, "\t", $altE-$altS+1, "\n";
					}
					
				}
				#New INV
				($invS, $invE) = ($info[0], $info[1]);
				($altC, $altS, $altE) = ($info[14], $info[2], $info[3]);
				$previousChr = $info[13];
				$previous = 1;
			}
			elsif($previous){
				#Continue INV
				($invE, $altE) = ($info[1], $info[3]);
			}
			else{
				#New INV
				$previous = 1;
				$previousChr = $info[13];
				($invS, $invE) = ($info[0], $info[1]);
				($altC, $altS, $altE) = ($info[14], $info[2], $info[3]);
			}
		}
		else{
			if($previous){
				#End INV
				$invLength += $invE - $invS + 1;
				my ($previousDiff, $postDiff) = getDiffs($previousChr, $invS-1, $invE+1);
				if($previousDiff =~ /[RQ]N/ || $postDiff =~ /[RQ]N/){
					$invGap++;
					$type = 'InversionGap';
				}
				else{
					$inv++;
					$type = 'Inversion';
					print "$type\t$previousChr\t$invS\t$invE\t$altC\t$altS\t$altE\t", $invE-$invS+1, "\t", $altE-$altS+1, "\n";
				}				
				($invS, $invE) = ();
			}
			$previous = 0;
			$previousChr = $info[13];
		}
	}
	close IN;
	return ($inv, $invGap, $invLength);
}

sub printTranslocation{
	open IN, $_coord_file;
	my ($transC, $transS, $transE, $trans, $transLength, $transGap, $type);
	my ($altC, $altS, $altE);
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		next if($l =~ /RN/);
		next if($l =~ /unanchor/);
		if($nameTable{$info[13]} ne $nameTable{$info[14]}){
		#if($info[13] ne $info[14]){
			if($info[13] eq $transC && $info[14] eq $altC){
				#Continue Trans
				($transE, $altE) = ($info[1], $info[3]);
			}
			else{
				#Translocation to Translocation
				if($altC){
					#End current translocation
					$transLength += $transE - $transS + 1;
					my ($previousDiff, $postDiff) = getDiffs($transC, $transS-1, $transE+1);
					if(!$previousDiff || !$postDiff){
						#print STDERR "C = $transC, S = $transS, E = $transE\n";
						#print STDERR "Previous = $previousDiff, Post = $postDiff\n";
					}
					
					if($previousDiff =~ /[QR]N/ || $postDiff =~ /[QR]N/){
						$transGap++;
						$type = 'TranslocationGap';
					}
					else{
						$trans++;
						$type = 'Translocation';
						print "$type\t$transC\t$transS\t$transE\t$altC\t$altS\t$altE\t", $transE-$transS+1, "\t", $altE-$altS+1, "\n";	
					}
				}
				#New Trans
				($transC, $transS, $altC, $altS, $transE, $altE) = ($info[13], $info[0], $info[14], $info[2], $info[1], $info[3]);
			}
		}
		else{
			if($altC){
				#End Trans
				$transLength += $transE - $transS + 1;
				my ($previousDiff, $postDiff) = getDiffs($transC, $transS-1, $transE+1);
				
				if($previousDiff =~ /[QR]N/ || $postDiff =~ /[QR]N/){
					$transGap++;
					$type = 'TranslocationGap';
				}
				else{
					$trans++;
					$type = 'Translocation';
					print "$type\t$transC\t$transS\t$transE\t$altC\t$altS\t$altE\t", $transE-$transS+1, "\t", $altE-$altS+1, "\n";	
				}				
				($transC, $transS, $transE, $altC, $altS, $altE) = ();
			}
		}
	}
	close IN;
	return ($trans, $transGap, $transLength);
}

sub count_SEQ{
	system("grep SEQ $_input > $_input.SEQ.txt");
	open IN, "$_input.SEQ.txt";
	my $count;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		my $chr = $info[0];
		$count++ if($info[5] eq $chr);
	}
	print "$_input: $count SEQ events\n";
}


sub GetSpecificSequences{
	open IN, $_input or die "Cannot open input file $_input\n";
	open TEMP, "+>$_input.temp.list";
	while(my $l = <IN>){
		chomp $l;
		my ($refCoord, $queryCoord) = split /\t/, $l;
		my ($qc, $qs, $qe) = $queryCoord =~ /(.*):(\d+)-(\d+)/;
		next if($qe - $qs + 1 < 100);
		print TEMP $queryCoord, "\n";
	}
	system("perl ~/scripts/fasta_tools.pl -T 4 -i $_fasta_file -j $_input.temp.list -o $_input.fasta");
	#system("rm $_input.temp.list");
}

sub tally_pav{
	die "Currently only supports reference based coordniates\n" if($_type eq 'Q' || $_type eq 'Query');
	open IN, $_list or die "Cannot open file list $_list\n";
	my @files;
	my $i = 0;
	while(my $l = <IN>){
		chomp $l;
		my ($f1, $f2) = split /\t/, $l;
		push @files, $f1;
		open FI, $f1;
		print STDERR "Tallying $f1, populating $f2\n";
		populate_coords($f2);
		$i++;
		while(my $line = <FI>){
			chomp $line;
			my @info = split /\t/, $line;
			next if($info[1] ne 'GAP');
			next if($info[4] < 100);
			if($_type eq 'Q' || $_type eq 'Query'){
				@info = convert_coordinate($line);
			}
			next unless @info;
			unless(update_overlap($info[0], $info[2], $info[3], $info[4], $i)){
				$PAV{$info[0]}{$info[2]}[0] = $info[3];
				$PAV{$info[0]}{$info[2]}[$i] = "$info[0]:$info[2]-$info[3]";
			}
		}
	}
	#Final check for internal overlap
	while(check_overlap()){
		print STDERR "Cleaning up overlaps in hash\n";
	}
	print STDERR "No more overlaps found\n";
	
	#Print Table
	print OUT "Coord\t", join("\t", @files), "\t";
	print OUT "Number of samples\tPAV size\n";
	foreach my $c(keys %PAV){
		#my @ary;
		foreach my $s(keys %{$PAV{$c}}){
			my @ary;
			my $count; 
			push @ary, "$c:$s-$PAV{$c}{$s}[0]";
			#print OUT "$c:$s-$PAV{$c}{$s}[0]\t";
			my $size = $PAV{$c}{$s}[0] - $s + 1;
			for(my $x = 1; $x <= @files; $x++){
				#print OUT $PAV{$c}{$s}[$x], "\t";
				push @ary, "$PAV{$c}{$s}[$x]";
				$count++ if $PAV{$c}{$s}[$x];
			}
			#print OUT $count, "\t", $size, "\n";
			push @ary, $count;
			push @ary, $size;
			print OUT join("\t", @ary), "\n";
		}
	}
	if($_get_seq){
		system("cut -f 1 $_output > $_output.temp.list");
		system("perl ~/scripts/fasta_tools.pl -T 4 -i $_fasta_file -j $_output.temp.list -o $_output.fasta");
		system("rm $_output.temp.list");
	}
}

sub check_overlap{
	foreach my $c(sort keys %PAV){
		foreach my $s(sort {$a<=>$b} keys %{$PAV{$c}}){
			my $e = $PAV{$c}{$s}[0];
			foreach my $start(sort {$a<=>$b} keys %{$PAV{$c}}){
				next if($s == $start);
				my $end = $PAV{$c}{$start}[0];
				if(get_max($start, $s) <= get_min($end, $e)){
					#Has overlap
					print STDERR "Invernal overlap: $c:$s-$e, $c:$start-$end\n";
					my ($new_start, $new_end) = (get_min($s, $start), get_max($end, $e));
					my @new_ary = merge_array(\@{$PAV{$c}{$s}}, \@{$PAV{$c}{$start}});
					$new_ary[0] = $new_end;
					#if($new_ary[13] > 150000){
					#	print STDERR "Something Wrong\n";
					#	for(my $i = 0; $i<=@files; $i++){
					#		print STDERR "$i = $files[$i]\t", $PAV{$c}{$s}[$i], "\t", $PAV{$c}{$start}[$i], "\t", $new_ary[$i], "\n";
					#	}
					#	<>;
					#}
					delete $PAV{$c}{$s};
					delete $PAV{$c}{$start};
					$PAV{$c}{$new_start} = [@new_ary];
					return 1;
				}
			}
		}
	}
	return 0;
}

sub update_overlap{
	my ($c, $s, $e, $size, $i) = @_;
	my ($new_start, $new_end);
	my $has_overlap = 0;
	if($PAV{$c}{$s}){
		#S already exist in hash
		$has_overlap = 1;
		$PAV{$c}{$s}[0] = get_max($PAV{$c}{$s}[0], $e);
		$PAV{$c}{$s}[$i] = "$c:$s-$e";
	}
	else{
		foreach my $start(keys %{$PAV{$c}}){
			my $end = $PAV{$c}{$start}[0];
			if(get_max($start, $s) < get_min($end, $e)){
				#has overlap, get the larger coordinate, update hash
				$has_overlap = 1;
				#This is the existing array
				my @tally = @{$PAV{$c}{$start}};
				if($s < $start){
					$new_start = $s;
					if($PAV{$c}{$s}){
						#Both start and s are in the hash
						#Merge tally
						@tally = merge_array(\@tally, \@{$PAV{$c}{$new_start}});
					}
					#Remove start
					delete $PAV{$c}{$start};
					#update
					$tally[0] = get_max($end, $e);	
					$PAV{$c}{$new_start} = [@tally];
					$PAV{$c}{$new_start}[$i] = "$c:$s-$e";
				}
				else{
					#Only keep start
					#update existing hash
					$new_start = $start;
					$tally[0] = get_max($end, $e);
					$PAV{$c}{$new_start}[$i] = "$c:$s-$e";
				}
			}
		}
	}
	return $has_overlap;
}

sub merge_array{
	my ($a1, $a2) = @_;
	#print STDERR "a1 = $a1, a2, $a2\n";
	my @ary1 = @$a1;
	my @ary2 = @$a2;
	#my (@ary1, @ary2) = (@$a1, @$a2);
	#print STDERR "ary1 = ", join(",", @ary1), "\n";
	#print STDERR "ary2 = ", join(",", @ary2), "\n";
	for(my $i = 1; $i <= @files; $i++){
		if($_details eq 'size'){
			#print STDERR "$i = $files[$i]\t$ary1[$i]\t$ary2[$i]\n";
			$ary1[$i] = 0 if(!$ary1[$i]);
			$ary2[$i] = 0 if(!$ary2[$i]);
			$ary1[$i] = max($ary1[$i], $ary2[$i]);
		}
		else{
			$ary1[$i] .= $ary2[$i];
		}
	}
	return @ary1;
=cut
	if(@ary1 > @ary2){
		for(my $i = 0; $i < @ary2; $i++){
			$ary1[$i] = $ary2[$i];
		}
		return @ary1;
	}
	else{
		for(my $i = 0; $i < @ary1; $i++){
			$ary2[$i] = $ary1[$i];
		}
		return @ary2;
	}
=cut
}

sub convert_coordinate{
	my $diff = shift;
	#print STDERR "Diff = $diff\n";
	my @info = split /\t/, $diff;
	#Get next block
	return 0 unless($coord_hash{$info[0]}{$info[3]+1});
	my @next_block = @{$coord_hash{$info[0]}{$info[3]+1}};
	#Get next block reference coordinate
	my ($nrc, $nrs, $nre) = $next_block[1] =~ /(\S+):(\d+)-(\d+)/;
	#Get prev block
	my ($pbc, $pbs, $pbe) = $next_block[2] =~ /(\S+):(\d+)-(\d+)/;
	my @prev_block = @{$coord_hash{$pbc}{$pbs}};
	my ($prc, $prs, $pre) = $prev_block[1] =~ /(\S+):(\d+)-(\d+)/;
	#Convert coordinate
	$info[2] = $pre + 1;
	$info[3] = $nrs - 1;
	return @info;
}


sub count_GAP{
	open IN, "$_output.GAP.sorted.out";
	my $count = 0;
	while(<IN>){
		chomp;
		my @info = split /\t/;
		$count++ if($info[4] > 100);
	}
	print "$_input", "\t", $count, "\n";
}

sub sort_GAP{
	open IN, $_input;
	my @prev_block;
	my @post_block;
	open OUT2, "+>temp.out";
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		next unless($info[1] eq 'GAP');
		#Get post block
		@post_block = @{$coord_hash{$info[0]}{$info[3]+1}};
		#print STDERR "Post block = ", join(";", @post_block), "\n";
		my ($post_qc) = $post_block[1] =~ /(\S+):/;
		#Get prev block
		my ($prc, $prs, $pre) = $post_block[2] =~ /(\S+):(\d+)-(\d+)/;
		@prev_block = @{$coord_hash{$prc}{$prs}};
		#print STDERR "Prev block = ", join(";", @prev_block), "\n";
		my ($prev_qc) = $prev_block[1] =~ /(\S+):/;
		if($info[0] eq $post_qc && $prc eq $prev_qc){
			print OUT2 "$l\n";
		}
		#<>;
	}
	system("sort -nr -k 7 temp.out > $_output.GAP.sorted.out");
}

sub sort_SEQ_INV{
	my $type = shift;
	my $open_case;
	my $last_case;
	my $pre_SV;
	my $post_SV;
	my $case_length;
	open IN2, $_coord_file;
	open OUT2, "+>temp.out";
	while(my $l = <IN2>){
		chomp $l;
		my @info = split /\t/, $l;
		my ($rc, $qc, $rd, $qd, $rs, $re) = ($info[13], $info[14], $info[11], $info[12], $info[0], $info[1]);
		
		my ($R, $Q)=
		$type eq 'INV' ? ($rd, $qd) : 
		$type eq 'SEQ' ? ($rc, $qc) : 
		('', '');
		
		#Get previous block
		my ($prc, $prs, $pre) = $coord_hash{$rc}{$rs}[2] =~ /(\S+):(\d+)-(\d+)/;
		
		if($R ne $Q){
			#INV or SEQ Case
			my @last_info = split /\t/, $last_case;
			my @case_info = split /\t/, $open_case;
			if($open_case){
				if($last_info[13] eq $rc && $last_info[14] eq $qc){
					#Continue same case
					$last_case = $l;
					$post_SV = $SV_hash{$rc}{$re+1};
				}
				else{
					#Close case
					my @pre_SV = split /\t/, $pre_SV;
					my @post_SV = split /\t/, $post_SV;
					$case_length = $last_info[1] - $case_info[0] + 1 + $pre_SV[4] + $post_SV[4];
					#print STDERR "open case = $open_case\nlast case - $last_case\nE = $last_info[1]\tS = $case_info[0]\nPre SV = $pre_SV\nPost SV = $post_SV\nPre SV size = $pre_SV[4]\tPost SV size = $post_SV[4]\nLength = $case_length\n";
					print OUT2 $case_info[13], "\t", $case_info[14], "\t", $case_info[0], "\t", $last_info[1], "\t", $case_length, "\n";
					#<>;
					#New case
					#print STDERR "Curr block = $rc:$rs-$re\nPrev Block = $prc:$prs-$pre\n";
					$open_case = $l;
					$last_case = $l;
					$pre_SV = $SV_hash{$prc}{$pre+1};
					$post_SV = $SV_hash{$rc}{$re+1};
				}
			}
			else{
				#New case
				#print STDERR "Curr block = $rc:$rs-$re\nPrev Block = $prc:$prs-$pre\n";
				$open_case = $l;
				$last_case = $l;
				$pre_SV = $SV_hash{$prc}{$pre+1};
				$post_SV = $SV_hash{$rc}{$re+1};
			}
		}
		else{
			if($open_case){
				#Close case
				my @last_info = split /\t/, $last_case;
				my @case_info = split /\t/, $open_case;
				my @pre_SV = split /\t/, $pre_SV;
				my @post_SV = split /\t/, $post_SV;
				$case_length = $last_info[1] - $case_info[0] + 1 + $pre_SV[4] + $post_SV[4];
				#print STDERR "open case = $open_case\nlast case - $last_case\nE = $last_info[1]\tS = $case_info[0]\nPre SV = $pre_SV\nPost SV = $post_SV\nPre SV size = $pre_SV[4]\tPost SV size = $post_SV[4]\nLength = $case_length\n";				
				print OUT2 $case_info[13], "\t", $case_info[14], "\t", $case_info[0], "\t", $last_info[1], "\t", $case_length, "\n";
				#<>;
				#Reset
				$open_case = '';
				$last_case = '';
			}
		}
	}
	close OUT2;
	system("sort -nr -k 5 temp.out > $_output.$type.sorted.out");
}

sub populate_coords{
	my $coord_file = shift;
	open CO, $coord_file;
	my ($prc, $prs, $pre);
	while(<CO>){
		chomp;
		my @info = split /\t/;
		my ($rc, $rs, $re);
		my ($qc, $qs, $qe);
		if($_type eq 'R' || $_type eq 'Ref'){
			($rc, $rs, $re) = ($info[13], $info[0], $info[1]);
			($qc, $qs, $qe) = ($info[14], $info[2], $info[3]);
		}
		elsif($_type eq 'Q' || $_type eq 'Query'){
			($rc, $rs, $re) = ($info[14], $info[2], $info[3]);
			($qc, $qs, $qe) = ($info[13], $info[0], $info[1]);	
		}
		
		#if($coord_hash{$rc}{$rs}){
		#	print STDERR "WARNING: Overlapping start coordinate for $rc:$rs\n";
		#}
		#if prev-coord, the next of prev-coord is the current coord.
		#only if ref chr matches
		if($prc && $prs && $pre){
			if($prc eq $rc){
				push @{$coord_hash{$prc}{$prs}}, "$rc:$rs-$re";
			}
			$coord_hash{$rc}{$rs} = [$re,"$qc:$qs-$qe", "$prc:$prs-$pre"];
		}
		else{
			$coord_hash{$rc}{$rs} = [$re,"$qc:$qs-$qe", ""];
		}
		($prc, $prs, $pre) = ($rc, $rs, $re);
		#print STDERR "PRC = $prc, PRS = $prs\n";
		#print STDERR "Current = $rc:$rs-$re\nPrevious = $coord_hash{$rc}{$rs}[2]\nNext from previous= @{$coord_hash{$prc}{$prs}}\n";
		#<>;
	}
	close CO;
}


sub populate_sv{
	open IN, $_input;
	while(my $l = <IN>){
		chomp $l;
		my @info = split /\t/, $l;
		$SV_hash{$info[0]}{$info[2]} = $l;
	}
	close IN;
}

sub get_max{
	my ($a, $b) = @_;
	return $a > $b ? $a : $b;
}

sub get_min{
	my ($a, $b) = @_;
	return $a > $b ? $b : $a;
}

sub max{
	my ($a, $b) = @_;
	return $a > $b ? $a : $b;
}

sub min{
	my ($a, $b) = @_;
	return $a > $b ? $b : $a;
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