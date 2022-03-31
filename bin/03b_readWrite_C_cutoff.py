#!/usr/bin/env python
#%%
import pysam
import glob, re, sys, os

#%%
sample_name = sys.argv[1]
path = sys.argv[2]
Ccutoff_list = {1:0, 3:0, 5:0, 15:0}


# Define functions
"""Count occurances of C's in bam file, remove them, write clean bam file to new file. """
def count_C(samfile):
    name = re.sub("_genomeMap_sorted.bam",'',samfile)
    sys.stdout.write("processing:\t" + os.path.basename(name)+'\n')
    removed_read_count = 0
    bamfile = pysam.AlignmentFile(samfile,"rb", check_sq=False, threads=4)

    # Open cutoff files for writing
    for k in Ccutoff_list.keys():
        save_name = name+"_"+str(k)+"_Ccutoff_PE.bam"
        print(save_name+'\n')
        Ccutoff_list[k] = pysam.AlignmentFile(save_name, "wb", template=bamfile)
    # check each read for C content, write to new file if passed filter
    for read in bamfile.fetch():
        C_count = 'null'
        sequence = read.query_sequence.upper()
        try:
            #if read.get_tag('YR') == 'G2A':  # read pair copy
             #   pass
            if read.get_tag('YG') == 'C2T':  # Aligned to positive strand
                C_count = float(sequence.count('C'))
            elif read.get_tag('YG') == 'G2A':  # Aligned to negative strand
                C_count = float(sequence.count('G'))
            else:
                print("tag not found")
                pass
        except KeyError:
            print("tag not found")
            break

        if C_count == 'null':
            pass
        elif C_count <= max(Ccutoff_list):
            for k,v in Ccutoff_list.items():
                if C_count <= k:
                    v.write(read)
        else:
            removed_read_count += 1
    bamfile.close()
    for k,v in Ccutoff_list.items():
        v.close()
    sys.stdout.write("\n" + " reads removed: " + str(removed_read_count))
    return removed_read_count


#%%
bam_file = path + "/working_data/" + sample_name + "/" + sample_name + "_genomeMap_sorted.bam"
print(bam_file)
count_C(bam_file)
