# MumSV

Collection of perl scripts to identify SVs from MUMer. 

Dependencies:
* Bio::SeqIO (for fasta_gap_info)


A typical workflow is as follows:

## Obtaining gap locations
`perl fasta_gap_info.pl -i <genome fasta file> > gapInfo.output.file`

## Annotate Gap to MUMer diff and coord output files
`perl MumSV.pl -T AnnotateGap -i <diff file> -j <coord file> -g <reference genome gapInfo> -G <query genome gapInfo> -c <Chromosome name correspondance table for refernece> -C <Chromosome name correspondance table for query>`

## NonGapAlignments
`perl MumSV.pl -T NonGapAlignments -i <Gap annotated diff file> -j <Gap annotated coord file>`
  
## Determine SV events
`perl MumSV.pl -T CountEvents -i <nonGap diff files> -j <nonGap coord files> -c <Chromosome name correspondance table for refernece> -C <Chromosome name correspondance table for query> `
