import itertools
from sys import argv, exit
from itertools import zip_longest

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


chunk_size=4


def trimandpaste(pathToFastaFile, output):
    #filename = pathToFastaFile.split('/')[-1]
    output = open(output,"w")
    with open(pathToFastaFile) as f:
        for lines in grouper(f, chunk_size, ""): #for every chunk_sized chunk
            header = lines[0]
            seq = lines[1]
            sep = lines[2]
            qual = lines[3]
            trimmed_seq = seq[2:-11]+seq[-6:-1]+"\n" # fooprint + barcode
            UMI = seq[0:2]+seq[-11:-6] #7nt in total 
            split_header = header.split(" ")
            new_header = split_header[0]+"_"+UMI+" "+split_header[1]
            if qual[-1:] == "\n":
                new_qual = qual[2:-11]+qual[-6:-1]+"\n"
            else:
                new_qual = qual[2:-10]+qual[-6:-1]
            output.write(new_header)
            output.write(trimmed_seq) 
            output.write(sep) 
            output.write(new_qual)

    output.close() 

def main():
    if len(argv) != 3: 
        exit("Usage: 2 arguments required\n1: Path to fasta file \n2: name of output file")

    # Get paths
    pathToFastaFile = argv[1]
    output = argv[2]
        
    trimandpaste(pathToFastaFile, output)

if __name__ == "__main__":
    main()
