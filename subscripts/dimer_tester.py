#!/usr/bin/python3

import sys, itertools, math, os
from pathlib import Path
from Bio import SeqIO

str_program = Path(__file__).name
dir_pd = Path(__file__).parent.parent
str_version = "1.3"

print("[dt] Starting honzo's {0} v{1} for dG/dTm estimation".format(str_program, str_version))
"""
HAHN'S NOTES
v1  - program now creates matrices for both delta-g and melting temperature comparisons
    - fixed extra indent where flt_tm_corrected is assigned only when priemrs are non-palindrome

"""

# default params
help_menu = \
"""
USAGE: python {p} <filename>

""".format(p=str_program)

# ======================================[ VARIABLE DELCARATIONS ]=====================================

global dic_ambs, complements, kmer_attr
dic_ambs = {"M":["A","C"], "R":["A","G"], "W":["A","T"], "S":["C","G"], \
    "Y":["C","T"], "K":["G","T"], "V":["A","C","G"], "H":["A","C","T"], \
    "D":["A","G","T"], "B":["C","G","T"], "N":["A","C","G","T"]}
complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
    'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K',
    'B':'V', 'V':'B', 'D':'H', 'H':'D', 'N':'N', "-":"N"}


# Getting User Arguments
if ("-h" or "--help") in sys.argv:
    print(help_menu); quit()
else:
    pth_input = Path(sys.argv[-1])

dir_work = pth_input.parent
pth_outdg = dir_work / (pth_input.stem + "_dg.csv")
pth_outtm = dir_work / (pth_input.stem + "_tm.csv")
pth_outsum = dir_work / (pth_input.stem + "_sm.csv")

print("[dt] Input file              : {0}".format(pth_input))
print("[dt] Output file (dg)        : {0}".format(pth_outdg))
print("[dt] Output file (tm)        : {0}".format(pth_outtm))
print("[dt] Output file (summary)   : {0}\n".format(pth_outsum))

# -------------------------------------[ Getting Pre-requisites ]-------------------------------------

class Kmer:
    def __init__(self, kmer, dg=-1):
        self.kmer = kmer
        self.dg = dg

    def set_tm(self, dh, ds, dg):
        self.tm_dh = dh
        self.tm_ds = ds
        self.tm_dg = dg

pth_kmerdg = dir_pd / "docs" / "experiment_2mers.csv"
pth_kmertm = dir_pd / "docs" / "santalucia-1996.csv"
kmer_attr = {}

print("[dt] Getting kmer delta-G profile from file: {0}".format(pth_kmerdg))
with open(pth_kmerdg, "rt") as rf_kmerdg:
    rf_kmerdg.readline()
    for line in rf_kmerdg:
        lst_index = line.strip().split(",")
        kmer_attr[lst_index[0]] = Kmer(lst_index[0], float(lst_index[1]))

# getting kmer information from docs (previous publication data)
print("[dt] Getting kmer Tm profile from file: {0}".format(pth_kmertm))
with open(pth_kmertm, "rt") as rf_kmertm:
    rf_kmertm.readline() # skip header line
    for line in rf_kmertm:
        lst_index = line.strip().split(",")
        kmer_attr[lst_index[0]].set_tm(dh=float(lst_index[1]), ds=float(lst_index[2]), dg=float(lst_index[3]))


# ============================================[ FUNCTIONS ]===========================================

def revc(seq):
    s = seq.upper().replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c")
    return s[::-1].upper()

def count_kmers(seq, k=2):
    dic_counts = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in dic_counts: dic_counts[kmer] += 1
        else: dic_counts[kmer] = 1

    return dic_counts

def get_gc(seq):
    # GC content should be the same for both sense and antisense
    seq = seq.upper()
    gc = (seq.count('C') + seq.count('G')) / len(seq) * 100
    return round(gc, 1)

def get_tm(seq, pconc=0.2e-6):
    # Using Nearest-Neighbour method from Santalucia-1996 to estimate melting temperature
    
    kmer_profile = count_kmers(seq)

    # initialization energies
    if seq.find("G") >= 0 or seq.find("C") >= 0:
        flt_dh = 0; flt_ds = -5.9; flt_dg = 1.82
    else:
        flt_dh = 0; flt_ds = -9.0; flt_dg = 2.8

    for kmer, count in kmer_profile.items():
        flt_dh += kmer_attr[kmer].tm_dh * count
        flt_ds += kmer_attr[kmer].tm_ds * count
        flt_dg += kmer_attr[kmer].tm_dg * count

    # convert from kcal/mol to cal/mol (eu)
    flt_dh *= 1000 
    if seq == revc(seq): # self-complementary strand
        flt_ds -= 1.4
        flt_dg += 0.4
        flt_tm = flt_dh / (flt_ds + 1.987 * math.log(pconc)) - 273.15
    else:
        flt_tm = flt_dh / (flt_ds + 1.987 * math.log(pconc/4)) - 273.15 
        # 0.797x + 6.8044
        #flt_tm_corrected = round(0.9834 * flt_tm - 17.36, 2)
    flt_tm_corrected = round(0.797 * flt_tm + 6.8044, 2)
    return flt_tm_corrected
    
def get_combinations(str_fullseq):
    """ Gets a sequence containing ambiguous bases, returns a list of sequences
    with all the possible ACTG combinations of that sequence """
    str_fullseq = str(str_fullseq)

    # we create a dictionary to figure out where all these ambiguous nucleotides are
    dic_amblocs = {}
    for i, char in enumerate(str_fullseq):
        if char in dic_ambs.keys():
            dic_amblocs[i] = char

    j = 0
    lst_segments = []
    #print("fullseq", str_fullseq)
    for i in dic_amblocs.keys():
        str_subseq = str_fullseq[j:i+1]
        #print(i, str_subseq)
        lst_subseqs = [str_subseq[:-1]+nuc for nuc in dic_ambs[str_subseq[-1]]]
        lst_segments.append(lst_subseqs)
        j = i+1

    str_endseq = str_fullseq[j:]
    lst_seqsegments = [seq for seq in itertools.product(*lst_segments)]
    lst_seqs = ["".join(seqsegment)+str_endseq for seqsegment in lst_seqsegments]
    return lst_seqs

def get_alignments(primer1, primer2):
    """
    Using the two primers, this function will shift base-by-base to model each
    possible configuration in the form of a string. Having generated the strings,
    est_binding() is called to estimate the strength of the attachment as well as
    other information returned with a list. 
    """

    def estimate_binding(sub_primer1, sub_primer2):
        """
        This function finds the strongest binding location from a given configuration
        passed through the function by sub_primer1 and sub_primer2. The result is an
        alignment string consisting of the following characters:
        ' ' not a complements between primers
        ':' complementary bases, but not part of the strongest connection, and
            therefore not used to calculate the binding strength
        '|' consecutive complementary bases used to calculate binding strength

        [input]                               [output]
        5'----GACTTACGTATT-----------3'       5'----GACTTACGTATT-----------3'
                                                    ||||     ::           
        3'-ACACTGAGCAGCTA------------5'       3'-ACACTGAGCAGCTA------------5'

        returns: [delta-G estimate, longest run of complements bases, alignment string]
        """
        dic_result = {}
        flt_currentest = 0
        int_run = 0
        flt_est = 999; int_endpos = -1; int_maxrun = -1
        str_aln = "" # will change with iteration

        for n in range(len(sub_primer1)):
            # check if it's a complements, if so, check if A-T or G-C
            if sub_primer1[n] == complements[sub_primer2[n]]:
                str_aln += ":"
                if sub_primer1[n] == "A" or sub_primer1[n] == "T": 
                    flt_currentest -= 2.0
                    int_run += 1
                elif sub_primer1[n] == "C" or sub_primer1[n] == "G": 
                    flt_currentest -= 3.0
                    int_run += 1

                if flt_currentest < flt_est:
                    flt_est = flt_currentest
                    int_endpos = n
                    int_maxrun = int_run
            else:
                str_aln += " "
                flt_currentest = 0
                int_run = 0

        # str_newaln will change the strongest continuous run of complementary
        # nucleotides (similar to Thermofisher's algorithm)
        str_newaln = str_aln[:int_endpos-int_maxrun+1] + "|"*int_maxrun + str_aln[int_endpos+1:]
        str_printout = "5'-"+sub_primer1+"-3'\n" \
            +"   "+str_newaln \
            +"\n3'-"+sub_primer2+"-5'\n"

        # get profile of 2-mers from primer1
        # the sub-sequence we're interested in should be at the same location as the pipes
        start = str_newaln.find("|")
        end = str_newaln.rfind("|")
        str_subseq1 = sub_primer1[start:end+1]
        dic_kmerprofile = count_kmers(str_subseq1)

        flt_estimate = 0
        for key in dic_kmerprofile:
            # estimate delta-G based on aligning sub-sequence: 10.1073/PNAS.83.11.3746
            flt_estimate -= kmer_attr[key].dg * dic_kmerprofile[key]
        
        # generated a regression line on excel, got the formula with R^2 = 0.9988
        flt_estimate = 1.0089 * flt_estimate + 0.0137
        dic_result["dg"] = round(flt_estimate, 2)
        dic_result["maxrun"] = int_maxrun
        dic_result["printout"]  = str_printout

        return dic_result

    # get_alignments function()
    dic_lowestresult = {"dg": 999}
    
    if len(primer1) > len(primer2):
        int_maxoffset = min(len(primer1), len(primer2)) - 1
    else:
        int_maxoffset = max(len(primer1), len(primer2)) - 1

    p1 = 0; p2 = 0
    # keey primer1 stationary, shift primer2 to the left
    sub_primer1 = primer1 + "-"*(int_maxoffset)
    sub_primer2 = "-"*(len(primer1)-1) + primer2 # not int_maxoffset
    while p2 < len(primer1):
        # estimate binding strength of this configuration
        # lst_result = [delta-G estimate, length of longest run, alignment string]
        dic_result = estimate_binding(sub_primer1, sub_primer2)
        if dic_result["dg"] < dic_lowestresult["dg"]:
            dic_lowestresult = dic_result.copy()

        sub_primer2 = sub_primer2[1:]+"-"
        p2 += 1

    # keep primer2 stationary, shift primer1 to the right
    sub_primer1 = "-" + primer1 + "-"*(int_maxoffset-1)
    sub_primer2 = primer2 + "-"*(len(primer1)-1)
    while p1 < len(primer2) - 1:
        #print(sub_primer1+"\n"+sub_primer2+"\n")
        dic_result = estimate_binding(sub_primer1, sub_primer2)
        if dic_result["dg"] < dic_lowestresult["dg"]:
            dic_lowestresult = dic_result.copy()
        sub_primer1 = "-"+sub_primer1[:-1]
        p1 += 1

    return dic_lowestresult

# =========================================[ START WORKFLOW ]=========================================

# get primer sequences from file, expand them and store in dictionaries
primer_seqs = {}
primer_tms = {}

print("[dt] Reading primer/IC sequences")
with open(pth_input) as handle:
    lst_records = list(SeqIO.parse(handle, "fasta"))

for record in lst_records:
    lst_combos = get_combinations(record.seq.upper())
    if len(lst_combos) <= 1:
        primer_seqs[record.id] = lst_combos[0]
        primer_tms[record.id] = get_tm(str(record.seq).upper())
    else:
        for i, str_expanded in enumerate(lst_combos, start=1):
            primer_seqs[record.id+"~"+str(i)] = str_expanded
            primer_tms[record.id+"~"+str(i)] = get_tm(str_expanded)

lst_pnames = list(primer_seqs.keys())  
    
print("[dt] Writing combination metrics to output files.")

with open(pth_outdg, "wt") as wf_outdg, open(pth_outtm, "wt") as wf_outtm, open(pth_outsum, "wt") as wf_outsum:
    # write header to files
    wf_outsum.write("name,sequence,tm,dg\n")
    wf_outdg.write("name,"+",".join(lst_pnames)+"\n")
    wf_outtm.write("name,"+",".join(lst_pnames)+"\n")

    for pname1, pseq1 in primer_seqs.items():
        # writing primer details to file
        flt_gc = get_gc(primer_seqs[pname1])
        flt_tm = primer_tms[pname1]
        aln_dg = get_alignments(pseq1, pseq1[::-1])
        flt_dg = aln_dg["dg"]
        wf_outsum.write("{n},{s},{t},{d}\n".format(n=pname1, s=pseq1, t=flt_tm, d=flt_dg))

        wf_outdg.write("{0}".format(pname1))
        wf_outtm.write("{0}".format(pname1))
        for pname2, pseq2 in primer_seqs.items():
            # all alignments should be 5>3 \n 3>5
            if pname1 == pname2:
                aln_result = aln_dg
            else:
                aln_result = get_alignments(pseq1, pseq2[::-1])

            wf_outdg.write(",{0}".format(aln_result["dg"]))
            flt_tmdiff = round(abs(primer_tms[pname2] - primer_tms[pname1]), 1)
            wf_outtm.write(",{0}".format(flt_tmdiff))

        wf_outdg.write("\n")
        wf_outtm.write("\n")