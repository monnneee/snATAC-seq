
import argparse
import pysam
import time, os
import gzip
from collections import defaultdict

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "input args")
    parser.add_argument("-fq","--fq",dest = "fastq",default = "")
    parser.add_argument("-b","--barcodewhitelist",dest = "barcode",default = "")
    parser.add_argument("-O", "--output", dest = "output",default = "barcode_correct.txt")

    return parser.parse_args()

def Mismatch(seq): #change the every base of seq to A&T&C&G (if seq length=20,there will be 80 new seq)
    base_list = ["A", "T", "C", "G"]
    seq_mut_list = []
    for i in range(len(seq)):
        seq_mut = [seq[:i] + base + seq[i+1:] for base in base_list]
        seq_mut = list(set(seq_mut)) #rm duplicate
        seq_mut_list = seq_mut_list + seq_mut
    return(seq_mut_list)

def GenerateMismatchDict(whitelist): #return all of barcode from whitelist to get (barcode_list); return a dict (barcode_dict) whose keys are all the mismatch(barcode) and values are all the original barcode.
    barcode_dict = defaultdict(set) #the dict class id set type
    filein = gzip.open( whitelist, "rt" )
    barcode_list = filein.readlines() 
    filein.close()
    barcode_list = [bc.strip() for bc in barcode_list]
    for line in barcode_list:
        barcode = line.strip().split("\t")[0]
        barcode_mut_list = Mismatch(barcode) 
        for seq in barcode_mut_list:
            if not barcode_dict[seq]:
                barcode_dict[seq].add(barcode)  
            else:
                continue
    for bc in barcode_list: #change the value of barcode_dict for the key with original barcode
        barcode_dict[bc] = {bc} 
    return(barcode_dict, barcode_list)

def main():

    parser = CommandLineParser()
    barcode_fq = parser.fastq
    barcode_lib = parser.barcode
    output = parser.output
    barcode_correct_file = output

    start_time = time.time()
    print("read whitelist",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    barcode_lib_dict, barcode_lib_list = GenerateMismatchDict(barcode_lib)
    end_time = time.time()
    print("read End", time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))

    start_time = time.time()
    print("correct ",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
    barcode_fq_in = pysam.FastxFile(barcode_fq)
    barcode_correct_file_out = open(barcode_correct_file, "w")
    for reads in barcode_fq_in:
        barcode = reads.name.split(":")[0][:28] #cut top 28 base of read name
        if barcode in barcode_lib_dict:
            for bc_value in list(barcode_lib_dict[barcode]):
                outstr = barcode + "\t" + "CB" + "\t" + bc_value + "\n" #if top 28bp of read name in the value list of barcode_dict(mismatch_barcode:barcode),print top28_readName and the matched value
                barcode_correct_file_out.write(outstr)
        else:
            continue
    barcode_correct_file_out.close()
    end_time = time.time()
    print("correct End", time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))

if __name__ == "__main__":
    main()


