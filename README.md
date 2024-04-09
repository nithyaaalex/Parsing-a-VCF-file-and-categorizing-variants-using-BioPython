Using BioPython to parse a GFF, VCF and a genome fasta file to categorize SNPs as non-coding, synonymous and non-synonymous

Requirements(libraries):
logging, argparse, os, gffutils, vcf, math, Bio.Seq, matplotlib.

The script must be run through a terminal or command prompt as it accepts 3 arguments using a command line interface, and all three are mandatory. 
It requires a vcf file to be provided using the --vcf named argument, a gff file using --gff and a fasta file using the --fasta argument. 

An example of how the script can be run is: 
python3 2933044.py --vcf assessmentData.vcf.gz --gffPlasmoDB-54_Pfalciparum3D7.gff --fasta PlasmoDB-54_Pfalciparum3D7_Genome.fasta


What the code does:
1. For each SNP in the VCF file, it checks if it lies inside a coding region. This is done by extracting all CDS that lies nearby and checking if the SNP is completely enveloped in it.
2. If it doesn't lie in one, it marks it as non-coding.
3. If it is, it translates the protein that is made by that gene originally and the new protein that occurs due to the SNP.
4. The two proteins are compared and if they are identical, it is marked as synonymous.
5. If not, it is marked as non-synonymous.
6. Outputs a file with these columns for each SNP in the VCF Chromsome,Position,Ref,Alt,Type(non-coding, synonymous or non-synonymous),Transcript,Ref AA,Alt AA
