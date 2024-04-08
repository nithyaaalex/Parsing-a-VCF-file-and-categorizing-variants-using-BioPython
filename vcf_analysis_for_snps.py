#import statements
import logging
import argparse
import os
import gffutils 
import vcf
import math
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from matplotlib import pyplot as plt



#setting up logger
mylogger = logging.getLogger()
mylogger.setLevel(logging.INFO) 
filehandle_log = logging.FileHandler('2933044_log_file.txt', mode = 'a') #choosing filehandler and giving the name of the log file
filehandle_log.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
mylogger.addHandler(filehandle_log) #adding filehandler to the logger created



#setting up argparse
parser_for_input = argparse.ArgumentParser(prog= 'Analysing SNPs to determine the effect the mutation has', description = 'Using the genome and fasta file provided, \
        we determine if the SNP is present in a noncoding, synonymous or non-synonymous region', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser_for_input.add_argument('--vcf', required= True, help= 'The name of the vcf file with the SNPs')
parser_for_input.add_argument('--gff', required= True, help= 'The name of the gff file for the same organism')
parser_for_input.add_argument('--fasta', required= True, help= 'The corresponding fasta file for the gff file')
args = parser_for_input.parse_args()



#checking if input provided is accessible
error = False #flag to exit the program incase files aren't found/ to log which files are found sucessfully
files_to_check = [args.vcf, args.gff, args.fasta]
for file in files_to_check:
    try:
        open(file)
    except FileNotFoundError:
        mylogger.error(f"The file:{file} can't be accessed, please try entering another file name for this format.\n")
        error = True

if error:
    raise SystemExit(1)
else:
    mylogger.info(f'The input files used are listed below:')
    mylogger.info(f'The VCF file used was: {args.vcf}')
    mylogger.info(f'The GFF file used was: {args.gff}')
    mylogger.info(f'The fasta file used was: {args.fasta}\n')



#function definitions
def cds_coordinates(variant_location, child_with_snp,  parent): 
    """
    function to find the cds location of the SNP
    it requires the location of the SNP, the cds that has the SNP and the parent transcript of the cds that has the SNP
    """
    if parent.strand == "+": #if else to set strandcheck accordingly for line 52
        strandcheck = False
    elif parent.strand == "-":
        strandcheck = True

    total = 0 #temporary variable to add the cds till the matching cds is found
    #iterating through every child of the transcript and adding the length of every cds till the cds with the snp is found 
    for child in db.children(parent.id, featuretype='CDS', order_by="start", reverse = strandcheck): 
        if child.id == child_with_snp:
            region_found = child
            break
        else:
            total += child.end - child.start + 1 #length of each cds is calculated end-start+1 (1 to account for subtraction removing 2 bases from either side)
 
    if parent.strand == "+": #calculating cds start if it's the positive strand - the total of each cds so far+ location of SNP - start of the cds with the snp +1
        relative_to_cds_start = total + variant_location - region_found.start + 1
    elif parent.strand == "-": #if it's the negative strand -  the total of each cds so far + the end of the cds with SNP - the location of snp +1
        relative_to_cds_start = total + region_found.end - variant_location + 1
    return relative_to_cds_start #returning the calculated value (it differs for each strand as the position values given are from 5' to 3' but the reverse has to be calculated in the opp direction)

def protein_coordinates(cds_coordinates): #divides cds coordinates by 3 and rounds upwards for protein position
    return math.ceil(cds_coordinates/3)

def check_synonymous(cds_position, protein_position, alternate_base, parent_of_cds, fasta_file):
    """
    this function checks if the SNP causes a synonymous change
    it requires the cds location of the SNP, the protein location of the change caused, -
    - the DNA base changed, parent of the cds in which the SNP was found and finally the fasta file in which the full sequence is stored
    """
    #this conditional sets the strandcheck and the base that changes according to strand
    if parent_of_cds.strand == "+":
        strandcheck = False
        base_for_substitution = str(alternate_base[0])
    elif parent_of_cds.strand == "-":
        strandcheck = True
        base_for_substitution = str(Seq((str(alternate_base[0]))).complement()) #complements the DNA base that changes for the reverse strand 

    seq = ''
    #for loop that constructes the whole cds sequence according to the strand it's on    
    for child in db.children(parent_of_cds.id, featuretype='CDS', order_by="start", reverse = strandcheck):
        seq = seq + child.sequence(fasta = fasta_file, use_strand=True)

    seq = Seq(seq) #making it a seq object to access its attributes
    mutated_seq = MutableSeq(seq) #creating a mutable version 
    mutated_seq[cds_position-1] = base_for_substitution #mutating the sequence at cds location -1 (due to python indexing)
    mutatedprotein_seq = Seq(mutated_seq).translate() #translating the mutated sequence
    protein_seq = seq.translate() #translating the original sequence

    #checking if both proteins generated are valid
    if (protein_seq.startswith("M") == True) and (protein_seq.endswith("*") == True) and (protein_seq.count("*") == 1): 
        if (mutatedprotein_seq.startswith("M") == True) and (mutatedprotein_seq.endswith("*") == True) and (mutatedprotein_seq.count("*") == 1):
            if protein_seq == mutatedprotein_seq: #synonymous mutation if there is no change in protein
                return (True, protein_seq[protein_position-1], "NA") #returns a flag that is used to determine synonymous/non-synonymous, the ref amino acid, and NA for alt amino acid
            else: #non-synonymous change if the proteins do not match 
                return (False, protein_seq[protein_position-1], mutatedprotein_seq[protein_position-1]) #returns a flag that is used to determine synonymous/non-synonymous, the ref amino acid and the alt amino acid
        else:
            return Exception #this exception would be caused by a mutation that caused the protein sequence to have an early stop codon/changed the start codon
    else:
        return Exception #this exception would be caused by pseudogenes



#main
try: #for logging success or failure of program
    vcfReader_input = vcf.Reader(filename= args.vcf) #using a variable to read in the vcf file
    db_name = args.gff.replace('gff', 'db') #changing the suffix of the gff file given to db 

    if not os.path.isfile(db_name): #checking if the database already exists in the path, if not it creates one
        mylogger.info("Creating database now ...\n")
        db = gffutils.create_db(args.gff, dbfn = db_name, keep_order = True, merge_strategy='merge', sort_attribute_values=True)
    else:
        mylogger.info("Connecting to existing database now ... \n") #if it already exists it connects to it 
        db = gffutils.FeatureDB(db_name, keep_order=True)

    #global count variables for bar plot/to track in logging
    count_quality = 0
    count_non_coding = 0
    count_synonymous = 0 
    count_non_synonymous = 0
    transcript_id = "" #needed global variable to change

    #opening output table file for
    output = open("2933044_bcpy_table.tsv", "w") 
    output.write(f"Chrom\tPos\tRef\tAlt\tType\tTranscript\tRef AA\tAlt AA\n")

    
    for record in vcfReader_input: #iterating over every SNP record
        if record.QUAL>20: #only considering those with quality higher than 20
            if record.REF != record.ALT: #only considering those where the base has actually changed
                region_found_flag = False #using a flag to signify if the record enters the if on line 155 - if it doesn't it is deemed non-coding

                #iterating over all the cds sequences that are on the same chromosome as the SNP and start lies near the snp position
                for feature in db.region(seqid=record.CHROM, start= record.POS, featuretype='CDS'):
                    if feature.start<=record.POS and feature.end>=record.POS: # only if the feature is completely enveloping the snp position the cds is identified
                        region_found_flag = True
                        for parent in db.parents(feature.id, featuretype='mRNA'): #iterating through all the parent transcripts of the cds
                            transcript_id = parent.id #saving the transcript id for table
                            cds_location = cds_coordinates(record.POS, feature.id ,parent) #finding the cds location of snp
                            protein_location = protein_coordinates(cds_location) #finding protein location of snp
                            try: 
                                check_change, ref_aa, alt_aa = check_synonymous(cds_location, protein_location, record.ALT, parent, args.fasta) #checking what type of change
                            except Exception:
                                mylogger.warning(f"an invalid protein was produced and this SNP was skipped")
                                break
                            
                            #assigning the variable of seq type according to the flags set, (true is synonymous, false is non-synonymous) and counting the total as well
                            if check_change == True :
                                seq_type = 'Synonymous'
                                count_synonymous += 1
                                break
                            else:
                                seq_type = 'Non-Synonymous'
                                count_non_synonymous += 1
                                break

                if region_found_flag == False: 
                    seq_type = 'Non-Coding' #non coding for the ones where cds wasn't found at that position
                    count_non_coding += 1
                    transcript = protein_location = ref_aa = alt_aa = "NA"   #setting the other variables to NA

                #writing output to the file 
                output.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{seq_type}\t{transcript_id}\t{protein_location}\t{ref_aa}\t{alt_aa}\n")
        else:
            count_quality+=1


    output.close() #closing the output file correctly

    #logging the no. of records that were of low quality
    mylogger.info(f'The number of records which did not pass the quality check  (quality<20) were : {count_quality}\n') 

    #code for creating the bar plot
    #setting the category names and values to show
    categories = ['Non-Coding', 'Synonymous', 'Non-Synonymous']
    values = [count_non_coding, count_synonymous, count_non_synonymous]
    plt.bar(categories, values, color='skyblue') #actually plots 

    #for loop for labelling the bars with the value 
    for index, value in enumerate(values): #enumerate returns the index and the value of the list values
        plt.text(index, value + 1, str(value), ha='center', va='bottom') #index is the x coordinate of where it is placed, value+1 is the y coordinate, str(value) is the actual value displayed, and ha and va help center it according to the bar (ha is horizontal axis)

    #labelling the plot and giving it a title
    plt.xlabel('Type of Sequence')
    plt.ylabel('Count')
    plt.title('Bar Plot comparison of the type of sequence the variant is present in(only variants with quality higher than 20)')
    plt.gcf().set_size_inches(10, 6) #increasing size of immage to accomodate for long title
    plt.savefig('2933044plot.png') #saving the plot

    #for telling the user where the output files have been stored
    current_directory = os.getcwd() #using os to get the current working directory
    mylogger.info(f"The program executed successfully, the output files have been stored in {current_directory}, they are called 2933044plot.png and 2933044_log_file.txt respectively.\n\n")

except Exception as e: #if any error was encountered, it will be reported in the log file
    mylogger.error(f"The program did not run correctly, the error/exception encountered was {e}\n\n")

logging.shutdown() #shutting logging down to close the logging file and exit the program properly