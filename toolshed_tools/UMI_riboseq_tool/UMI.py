import gzip
from mimetypes import guess_type
from functools import partial
from sys import argv, exit
import itertools
from itertools import zip_longest
import subprocess
from subprocess import call
import Bio
from Bio import SeqIO

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def UMI_processing(pathToFastaFile, output_path):
    
    output = open(output_path,"w")
    
    if is_gz_file(pathToFastaFile) == True: 
        print ('file is gzipped fastq')
        with gzip.open(pathToFastaFile, "rt") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                lines = record.format('fastq').split('\n') # all 4 lines 
                header = lines[0]
                if i % 100000 == 0:
                    print ('read number %s' % i)
                seq = lines[1]
                sep = lines[2]
                qual = lines[3]
                if (header.startswith('@')):
                    trimmed_seq = seq[2:-6] # fooprint + barcode
                    UMI = seq[0:2]+seq.rstrip()[-5:len(seq)] #7nt in total; 5'NN and last 3'NNNNN
                    split_header = header.split(" ")
                    new_header = split_header[0]+"_"+UMI+" "+split_header[1]
                    new_qual = qual[2:-6]
                    output.write(new_header+'\n')
                    output.write(trimmed_seq+'\n')
                    output.write(sep+'\n')
                    output.write(new_qual+'\n')   
                    
    else: 
        for record in SeqIO.parse(pathToFastaFile, 'fastq'):
            lines = record.format('fastq').split('\n') # list of each record: id, seq, '+', quality 
            header = lines[0]
            seq = lines[1]
            sep = lines[2]
            qual = lines[3]
            trimmed_seq = seq[2:-6] # fooprint + barcode
            UMI = seq[0:2]+seq.rstrip()[-5:len(seq)] #7nt in total; 5'NN and last 3'NNNNN  
            split_header = header.split(" ")
            new_header = split_header[0]+"_"+UMI+" "+split_header[1]
            new_qual = qual[2:-6]
            output.write(new_header+'\n')
            output.write(trimmed_seq+'\n') 
            output.write(sep+'\n') 
            output.write(new_qual+'\n')
            
    output.close()
            

    

def main():
    if len(argv) != 3:
        exit("Usage: 2 arguments required\n1: Path to fasta file \n2: name of output file")

    # Get paths
    pathToFastaFile = argv[1]
    output = argv[2]
    UMI_processing(pathToFastaFile, output)

if __name__ == "__main__":
    main()
