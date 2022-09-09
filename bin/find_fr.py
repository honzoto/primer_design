import sys, os, re, glob
from pathlib import Path
import logging



"""
PURPOSE: To look for strain-specific primer sequences using a common probe

Normally, after running filter_db, we would run generate_psets.py to find potential combinations of f/p/r primers.
However, an alternative idea is to look for common probes, and then to generate strain-specific primers.

To do this, we have to first run get_primers.py and filter_db,py on all the strains we are interested in.
The output is a list of filtered primer sequences for each strain.

First, we have to perform an algorithm similar to get_locations.py to determine where in the genome the probe

"""

str_program = Path(__file__).name
dir_program = Path("/mnt/tank/bench/scripts/primer_design")
str_version = "1.0.1"

help_menu = """
----------------------------------[ HELP MENU ]---------------------------------

    -h/--help       : shows this menu
    -p/--project    : <dir> name of the project folder
    -i/--input      : <path> name of the FASTA file with filtered probes
    -t/--targets    : <path> name of the metadata (csv) file containing info
                        on which FASTA genome is matched with which filtered
                        primer sequences FASTA

--------------------------------------------------------------------------------

"""

# ======================================[ VARIABLE DECLARATIONS ]=====================================

global complements
lst_amplimits = [50, 200]
flt_mindg = -12.0
lst_tmprobes = [-3.0, -8.0]
complements = {"A":"T", "T":"A", "C":"G", "G":"C", "-":"N"}

# get dimer profile from file
pth_kmerdg = Path("docs") / "experiment_2mers.csv"
pth_kmertm = Path("docs") / "santalucia-1996.csv"

# ------------------------------------[ Retrieving User Arguments ]-----------------------------------

if len(sys.argv) > 1:
    for i, arg in enumerate(sys.argv):
        if arg == ("-p" or "--project"):
            str_project = sys.argv[i+1]
        elif arg == ("-i" or "--input"):
            str_probefile = sys.argv[i+1]
        elif arg == ("-t" or "--targets"):
            str_targetfile = sys.argv[i+1]

        elif arg == ("-h" or "--help"):
            print(help_menu); quit()
else:
    print("[ERROR] no arguments specified.")
    print(help_menu); quit()


dir_project = Path("/mnt/tank/bench/projects") / str_project
pth_probe = dir_project / str_probefile
pth_targets = dir_project / str_targetfile
pth_log = dir_project / "PrimerDesign.log"

# Setting up logger
print("[ff] Setting up logger: {0}".format(pth_log))
#open(pth_log, "wt").close() # reset the log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
log_formatter = logging.Formatter("%(asctime)s-%(levelname)s: %(message)s", datefmt="[%Y-%m-%d %H:%M:%S]")
file_handler = logging.FileHandler(pth_log)
file_handler.setFormatter(log_formatter)
logger.addHandler(file_handler)


# ======================================[ CLASSES AND FUNCTIONS ]=====================================

# CLASSES
class Kmer:
    def __init__(self, kmer, dg=-1):
        self.kmer = kmer
        self.dg = dg

    def set_tm(self, dh, ds, dg):
        self.tm_dh = dh
        self.tm_ds = ds
        self.tm_dg = dg

# FUNCTIONS
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
        logger.warning("Cannot shorten primer name: {0}".format(expanded_name))
        return expanded_name

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
        passed through the function by sub_primer1 and sub_primer2. 
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

def ham_dist_ambs(seq1, seq2):
    ambs = {"A": ["A"], "C":["C"], "T":["T"], "G":["G"], "-":["-"],
        "M":["A","C"], "R":["A","G"], "W":["A","T"], "S":["C","G"],
        "Y":["C","T"], "K":["G","T"], "V":["A","C","G"], "H":["A","C","T"],
        "D":["A","G","T"], "B":["C","G","T"], "N":["A","C","G","T"]}
    ambs_rev = {"A":"A", "C":"C", "T":"T", "G":"G",
        "AC":"M", "AG":"R", "AT":"W", "CG":"S", "CT":"Y", "GT":"K",
        "ACG":"V", "ACT":"H", "AGT":"D", "CGT":"B", "ACGT":"N"}

    # arguments entered as (search, source)
    if len(seq1) != len(seq2):
        raise ValueError("Undefined")

    hd = 0
    for nuc1, nuc2 in zip(seq1, seq2):
        com1 = set(ambs[nuc1])
        com2 = set(ambs[nuc2])
        hd += len(com1^com2) / len(com1|com2)

    return round(hd, 2)

# ==========================================[ MAIN WORKFLOW ]=========================================

# header block
print("[gs] Starting {0} for primer attribute analysis".format(str_program))
str_log = "\n--- STARTING [{0}] FOR PRIMER ATTRIBUTE ANALYSIS ---\n".format(str_program)
str_log += "\tPrimer input file:                            : {0}\n".format(pth_probe.name)

logger.info(str_log+"\n")

# ------------------------------------[ Getting input information ]-----------------------------------


# getting kmer information from docs (previous publication data)
global kmer_attr, target_pairs
kmer_attr = {}
target_pairs = {}

print("[ff] Getting kmer Tm profile from file: {0}".format(pth_kmertm))
with open(pth_kmertm, "rt") as rf_kmertm:
    rf_kmertm.readline() # skip header line
    for line in rf_kmertm:
        lst_index = line.strip().split(",")
        kmer = Kmer(lst_index[0])
        kmer.tm_dh, kmer.tm_ds, kmer.tm_dg = float(lst_index[1]), float(lst_index[2]), float(lst_index[3])
        kmer_attr[lst_index[0]] = kmer

print("[ff] Getting kmer delta-G profile from file: {0}".format(pth_kmerdg))
with open(pth_kmerdg, "rt") as rf_kmerdg:
    rf_kmerdg.readline()
    for line in rf_kmerdg:
        lst_index = line.strip().split(",")
        kmer_attr[lst_index[0]].dg = float(lst_index[1])

print("[ff] Getting target pairs from file: {0}".format(pth_targets))
with open(pth_targets, "rt") as rf_targets:
    rf_targets.readline()
    for line in rf_targets:
        lst_index = line.strip().split(",")
        # primer sequences: reference consensus genome
        target_pairs[lst_index[0]] = lst_index[1]
        

# Get locations of the probes on each read
# get_consensus.py should have generated a consensus sequence from each of the target genomes
# the genomes should be paired up with a filtered primer file
for str_primerfile, str_referencefile in target_pairs.items():
    print("[ff] Processing primers from {0} with {1}...".format(str_primerfile, str_referencefile))
    



print("[ff] Program complete.")