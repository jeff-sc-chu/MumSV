# MumSV

Collection of perl scripts to identify SVs from MUMer. 

A typical workflow is as follows:
#Obtaining gap locations
perl fasta_gap_info -i <genome fasta file> > gapInfo.output.file

# Annotate Gap to MUMer diff and coord output files
perl ~/scripts/analyze_mummer.pl -T AnnotateGap -i <diff file> -j <coord file> -g <reference genome gapInfo> -G <query genome gapInfo> -c <Chromosome name correspondance table for refernece> -C <Chromosome name correspondance table for query>

# NonGapAlignments
perl ~/scripts/analyze_mummer.pl -T NonGapAlignments -i <Gap annotated diff file> -j <Gap annotated coord file>
  
# Count events
perl ~/scripts/analyze_mummer.pl -T CountEvents -i <nonGap diff files> -j <nonGap coord files> -c <Chromosome name correspondance table for refernece> -C <Chromosome name correspondance table for query> 
