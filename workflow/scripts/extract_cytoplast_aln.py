import pysam
#tab separated table, first column: sample name, second column: path to bam file
sample_table="samples.tsv"
#minimum sequencing depth
min_depth=5
#name of contigs we want to look at, e.g. for cpDNA
sample_contigs=["OY340914.1"]
#name of the output file
output_aln="cp_alignments.fasta"
#dictionary to hold results
aln_dict={}
with open(sample_table) as sample_files:
    for line in sample_files:
        sample_name=line.strip().split()[0]
        sample_bam=line.strip().split()[1]
        print(sample_name)
        #Use pysam to load alignment file
        bam_input=pysam.AlignmentFile(sample_bam, "rb")
        #make empty string to add output to
        aln_string=""
        #Do a loop, because we can have multiple contigs
        for sample_contig in sample_contigs:
            #do pileup
            bam_pileup=bam_input.pileup(sample_contig)
            record_depth_sum=0
            #iterate over each column
            for pileup_pos in bam_pileup:
                #test if we have minimum depth at position, otherwise just add N to output sequence
                if pileup_pos.n<=min_depth:
                    aln_string=aln_string+"N"
                else:
                    pos_count={"A":0, "C":0, "G":0, "T":0, "N":0}
                    #Here is the tricky part, I think the refskip part is more complex
                    for pileupread in pileup_pos.pileups:
                        if pileupread.is_del or pileupread.is_refskip:
                            pos_count["N"]+=1
                        else:
                            pos_count[pileupread.alignment.query_sequence[pileupread.query_position]]+=1

                    aln_string=aln_string+[w for w in sorted(pos_count, key=pos_count.get, reverse=True)][0]
        aln_dict.update({sample_name:aln_string})
#write results for all samples to one fasta file
with open(output_aln, "w") as output_file:
    for output_sample in aln_dict:
        output_file.write(">%s\n%s\n" % (output_sample, aln_dict[output_sample]))
