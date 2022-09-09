#!/home/sbsuser/miniconda3/bin/python3

import sys, glob, os, math, re
from pathlib import Path

version = "1.0"
program = Path(__file__).name
print("This is Honzo's {0} v{1} for probe-based primer finding".format(program, version))
"""
This program takes a probe sequence and looks for potential primer sequences from multiple genomes
that will fit the amplicon criteria


"""

# ======================================[ VARIABLE DECLARATIONS ]=====================================

"""
-------------------------------< HELP MENU >------------------------------------

    -i/--input          : <path> FASTA file with possible probe seqs
    -r/--primer         : <path> FASTA file with potential primer seqs

--------------------------------------------------------------------------------
"""

for i, arg in enumerate(sys.argv):
    if arg == ("-p" or "--probe"):
        str_project = sys.argv[i+1]
    elif arg == ("-i" or "--input"):
        str_input = sys.argv[i+1]
    elif arg[0] == "-":
        print("[WARNING] unrecognized argument {0}".format(arg))
    

dir_pd = Path("/mnt/tank/bench/scripts/primer_design")


global kmer_attr, complements
complements = {"A":"T", "T":"A", "C":"G", "G":"C", "-":"N"}


# ======================================[ CLASSES AND FUNCTIONS ]=====================================

class Kmer:
    def __init__(self, kmer, dg=-1):
        self.kmer = kmer
        self.dg = dg

    def set_tm(self, dh, ds, dg):
        self.tm_dh = dh
        self.tm_ds = ds
        self.tm_dg = dg

class Primer:
    """Each primer will have intrinsic characteristics stored in this object
    for record-keeping purposes, all primers sequences go 5>3 on the sense strand
    unless otherwise indicated through the .strand attribute"""

    def __init__(self, id: str, seq):
        #seq = seq.upper()
        self.id = id
        self.strand = "sense"
        self.id_com = shorten_name(self.id)
        self.seq = str(seq).upper().replace("U","T")
        self.len = len(self.seq)
        self.is_sensical = True # hz191
        self.revc = revc(self.seq)
        self.get_coords()

    def get_coords(self):
        for key, value in get_fields(self.id_com).items():
            setattr(self, key, value)
        
    def get_gc(self):
        # GC content should be the same for both sense and antisense
        gc = (self.seq.count('C') + self.seq.count('G')) / self.len * 100
        self.gc = round(gc, 1)
        return self.gc

    def get_tm(self, pconc=0.2e-6):
        # Using Nearest-Neighbour method from Santalucia-1996 to estimate melting temperature
        
        # Tm can be different for sense and antisense
        seqs_fr = {"tm": self.seq, "tm_revc": self.revc}
        for strand, seq in seqs_fr.items():
            kmer_profile = count_kmers(seq)
            if seq.find("G") >= 0 or seq.find("C") >= 0:
                flt_dh = 0; flt_ds = -5.9; flt_dg = 1.82
            else:
                flt_dh = 0; flt_ds = -9.0; flt_dg = 2.8

            for kmer, count in kmer_profile.items():
                flt_dh += kmer_attr[kmer].tm_dh * count
                flt_ds += kmer_attr[kmer].tm_ds * count
                flt_dg += kmer_attr[kmer].tm_dg * count
            flt_dh *= 1000 # convert from kcal/mol to cal/mol (eu)

            if self.seq == self.revc: # self-complementary strand
                flt_ds -= 1.4
                flt_dg += 0.4 # for symmetry correction
                flt_tm = flt_dh / (flt_ds + 1.987 * math.log(pconc)) - 273.15
            else:
                flt_tm = flt_dh / (flt_ds + 1.987 * math.log(pconc/4)) - 273.15

            # perform correction for salt concentrations (see tm_predictions.xlsx)
            flt_tm_corrected = round(0.797 * flt_tm + 6.8044, 2)
            setattr(self, strand, flt_tm_corrected)

        # hz249 - we set tm, tm_revc attributes using setattr function
        return self.tm, self.tm_revc

    def get_dg(self):
        results = get_alignments(self.seq, self.seq[::-1])
        self.dg = results["dg"]
        return self.dg

def revc(seq, rev=True):
    if rev: return "".join(complements.get(base, base) for base in reversed(seq))
    else: return "".join(complements.get(base, base) for base in seq)
        
def count_kmers(seq, k=2):
    dic_counts = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in dic_counts: dic_counts[kmer] += 1
        else: dic_counts[kmer] = 1
    return dic_counts

def get_fields(primer_id):
    result = {}
    m = re.search('^>?(.+)\|(\d+)\-(\d+)\|(\d+)\|(\d+)\~?(\d+)?', primer_id)
    result["gene"] = m.group(1)
    result["start"] = int(m.group(2))
    result["end"] = int(m.group(3))
    result["mm"] = int(m.group(4))
    try: result["variant"] = m.group(5)
    except: pass
    return result
    
def shorten_name(expanded_name):
    try: # see if there's an easier way to do this regex search
        m = re.search('^>?([^\|]*)\|([^\|]*)\|(\d+)\|(\d+)(\~\d*)?', expanded_name)
        return "{0}|{1}|{2}|{3}".format(m.group(1), m.group(2), m.group(3), m.group(4))
    except:
        print("Cannot shorten primer name: {0}".format(expanded_name))
        return expanded_name

# 1. For Delta-G relative estimations
def get_alignments(primer1, primer2):
    """
    Using the two primers, this function will shift base-by-base to model each
    possible configuration in the form of a string. Having generated the strings,
    est_binding() is called to estimate the strength of the attachment as well as
    other information returned with a list. """

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
    dic_lowestresult = {"dg":999}
    
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



# -----------------------------------[ Getting constants from lib ]-----------------------------------
pth_kmerdg = Path("docs") / "experiment_2mers.csv"
pth_kmertm = Path("docs") / "santalucia-1996.csv"

kmer_attr = {}
# getting kmer information from docs (previous publication data)
print("[gs] Getting kmer Tm profile from file: {0}".format(pth_kmertm))
with open(pth_kmertm, "rt") as rf_kmertm:
    rf_kmertm.readline() # skip header line
    for line in rf_kmertm:
        lst_index = line.strip().split(",")
        kmer = Kmer(lst_index[0])
        kmer.tm_dh, kmer.tm_ds, kmer.tm_dg = float(lst_index[1]), float(lst_index[2]), float(lst_index[3])
        kmer_attr[lst_index[0]] = kmer

print("[gs] Getting kmer delta-G profile from file: {0}".format(pth_kmerdg))
with open(pth_kmerdg, "rt") as rf_kmerdg:
    rf_kmerdg.readline()
    for line in rf_kmerdg:
        lst_index = line.strip().split(",")
        kmer_attr[lst_index[0]].dg = float(lst_index[1])


# ======================================[ WORKFLOW STARTS HERE ]======================================

all_probes = {}
all_primers = {}

