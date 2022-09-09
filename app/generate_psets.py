#!/usr/bin/python3

import sys, os, re, math
from pathlib import Path
import logging

from Bio import SeqIO

"""
PURPOSE: 
This program takes a list of possible primers saved in a FASTA file,
filters out poor primers based on GC content, 

HAHN's NOTES:
v1  - for a primer set to be considered valid, we need the binding strength
        for all ambiguous nucleotides - our test was made with 0 mismatches
    - added amplicon limit, tm-difference limits to filtering options
    - we have to check that gc content option has been selected before filtering
        otherwise it will cause an error - done
v2  - added step to check for heterodimers for internal controls, and remove them
        from the analysis
    - the alignment graphic now indicates 5'-3' and 3'-5' notation. Before this,
        we were aligning the actual sequences together which doesn't make sense
v3  - new delta-G estimation algorithm (from IDT study)
    - added check reference option to see if primer to eliminate is actually part of a reference
    - can now check for forward/probe/reverse all in one script, and verify dG for all three
    - storing primer combinations in Pcombo class, so details can be accessed later for write
    - changed primer1/primer3/primer2 notation to primer_f/primer_p/primer_r for ease of understanding
v4  - new structure for taqman primer sets; checking between sense/antisense for probe
    - detailed diagnostics for why primer combos don't make the cutoff (in dictionary)
    - completed first functional version of the pipeline for both TM and EP primer designs
    - hz191 - made check_reference part of Primer class, every primer will now have an is_sensical attr
        - if we don't have a reference, all primers are considered sensical, hence we reduce line numbers
    - created getattr function for Pcombo class to output mean dG and Tm to *_sets.csv output
    - psets now has individual primer set characteristics - as per michael's comments
        - program also writes pfails for primer sets that did not meet pset criteria (in same format)
    - created n_primers and n_sensicals to Pcombo class, so we know how many primers are in solution
v5  - when checking for Tm limits, we will only eliminate it if both the sense and antisense fails
        - this is because at the preliminary filtering stage, we don't know if a primer will be fwd/rev
    - hz627 added internal control analysis to get_alignments function
    - hz916 not writing alignment to file to save space
    - added option to write primers only to fasta, disregarding combinations

PROBLEMS
this script is not generating the appropriate amplicon sizes, for example:
aligned|95-118|24|0	55.5	-6.95	aligned|245-264|20|1	sense	58	-6.34	aligned|1533-1556|24|0	55.1	-3.62	3	118	-0.4	-5.16	2.7


"""


# ======================================[ VARIABLE DECLARATIONS ]=====================================

str_version = "1.5.3"
str_program = "generate_psets.py"
dir_pd = Path("/mnt/tank/bench/scripts/primer_design")
os.chdir(dir_pd)
print("[gs] Starting {p} (version {v}) for heterodimer estimation\n".format(p=str_program, v=str_version))

# we have to assign global names before declaring the variables
global complements, gene_seqs, kmer_attr
complements = {"A":"T", "T":"A", "C":"G", "G":"C", "-":"N"}

# get dimer profile from file
pth_kmerdg = Path("docs") / "experiment_2mers.csv"
pth_kmertm = Path("docs") / "santalucia-1996.csv"
dir_ic = Path("internal_controls")

global flt_mindg, flt_maxgcdiff, flt_maxtmdiff, lst_tmprobes

# Setting default parameters (some can be overwritten by sys.argv)
tf_primersonly = False

lst_runtype = ['dg', 'gc', 'tm', 'ic']
lst_gclimits = [30, 70] # [min, max] GC content for primer pool
lst_tmlimits = [55, 65] # [min, max] Tm estimate for primer pool
lst_tmprobes = [2.0, 6.0] # prb Tm must be higher than average of fwd/rev by [min, max]

flt_mindg = -10.0 # minimum delta G for primer combination
flt_maxgcdiff = 999 # %GC between fwd/rev must not exceed this
flt_maxtmdiff = 4.0 # Tm between fwd/rev must not exceed this
tf_checkic = False # Check alignments for interactions with internal controls

str_reference = "none" # if we discover a negative primer that's not part of the genome, does not matter
str_kit = "taqman"
lst_amplimits = [50, 180]

help_menu = \
"""
----------------------------------[ HELP MENU ]---------------------------------

    -h/--help       : shows this menu
    -p/--project    : <directory> name of the project directory
    -i/--input      : <file> name of the FASTA file containing the primer
                        sequences that are to be analyzed
    -k/--kit        : [default={k}] type of kit (taqman/endpoint)
    -r/--reference  : [default={r}] indicate reference genomes (fasta) to 
                        check if a negative result affects the binding
    --gc_limits     : [default={gc}] only keep primers with GC content
                        between and including lower and upper limits
    --tm_limits     : [default={tm}] only keep primers with basic Tm
                        values within the specified range
    --dg_limit      : [default={dg}] minimum delta-G value cutoff
    --primers_only  : [default={po}] stop program after performing filtration

--------------------------------------------------------------------------------
""".format(k=str_kit, gc=lst_gclimits, tm=lst_tmlimits, dg=flt_mindg, r=str_reference, po=tf_primersonly)

# -----------------------------------[ Getting arguments from user ]----------------------------------

if len(sys.argv) > 1:
    for i, arg in enumerate(sys.argv):
        if arg == ("-p" or "--project"):
            str_project = sys.argv[i+1]
        elif arg == ("-i" or "--input"):
            str_input = sys.argv[i+1]
        elif arg == ("-k" or "--kit"):
            str_kit = sys.argv[i+1]
        elif arg == ("-r" or "--reference"):
            str_reference = sys.argv[i+1]

        # filtering parameters
        elif arg == "--primers_only":
            tf_primersonly = True
        elif arg == "--gc_limits":
            lst_gclimits_ = sys.argv[i+1].split(",")
            lst_gclimits = [float(v.strip()) for v in lst_gclimits_]
        elif arg == "--tm_limits":
            lst_tmlimits_ = sys.argv[i+1].split(",")
            lst_tmlimits = [float(v.strip()) for v in lst_tmlimits_]
        elif arg == "--dg_limit":
            # if user enters positive dG, value is coerced to negative
            flt_mindg = abs(float(sys.argv[i+1])) * -1
        # generic flags
        elif arg == ("-h" or "--help"):
            print(help_menu)
            quit()
        elif arg[0] == "-":
            print("[warning] unrecognised flag:", arg)

else:
    print(help_menu)
    quit()

# -----------------------------[ Setting directories and user variables ]-----------------------------

dir_project = Path("/mnt/tank/bench/projects") / str_project
pth_input = dir_project / str_input

pth_log = dir_project / "PrimerDesign.log"
pth_out_primers = dir_project / (pth_input.stem + "_passing.fasta")
pth_out_aln = dir_project / pth_input.name.replace(".fasta", "_{0}_alns.txt".format(str_kit))
pth_out_psets = dir_project / pth_input.name.replace(".fasta", "_{0}_sets.csv".format(str_kit))
pth_out_fails = dir_project / pth_input.name.replace(".fasta", "_{0}_failsets.csv".format(str_kit))

if str_kit.lower() == "taqman": 
    lst_amplimits = [60, 150]
    flt_uppertm = lst_tmlimits[1] + sum(lst_tmprobes) / 2
    lst_tmlimits[1] = flt_uppertm

elif str_kit.lower() == "endpoint":
    lst_amplimits = [300, 400]

# setting up logger
print("[gs] Setting up logger: {0}".format(pth_log))
# open(pth_log, "wt").close() # reset the log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
log_formatter = logging.Formatter("%(asctime)s-%(levelname)s: %(message)s", datefmt="[%Y-%m-%d %H:%M:%S]")
file_handler = logging.FileHandler(pth_log)
file_handler.setFormatter(log_formatter)
logger.addHandler(file_handler)

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

    def check_reference(self):
        self.is_sensical = self.seq in "\t".join([str(r).upper() for r in gene_seqs.values()])
        return self.is_sensical
        
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
            #flt_tm_corrected = round(0.9834 * flt_tm - 17.36, 2)
            flt_tm_corrected = round(0.797 * flt_tm + 6.8044, 2)
            setattr(self, strand, flt_tm_corrected)

        # hz249 - we set tm, tm_revc attributes using setattr function
        return self.tm, self.tm_revc

    def get_dg(self):
        results = get_alignments(self.seq, self.seq[::-1])
        self.dg = results["dg"]
        return self.dg

class Pcombo:
    def __init__(self, id_f: str, id_r: str, primers_f: list, primers_r: list, kit="endpoint"):
        self.id_f = id_f
        self.id_r = id_r
        self.primers_f = primers_f
        self.primers_r = primers_r
        self.kit = kit

    def add_probe(self, id_p: str, primers_p: list, strand_p="sense"):
        self.id_p = id_p
        self.primers_p = primers_p
        self.strand_p = strand_p

    def get_ampsize(self) -> int:
        self.fields_f = get_fields(self.id_f)
        self.fields_r = get_fields(self.id_r)

        if self.kit == "endpoint": 
            # if the amplicon region is too large, the loop calling this fn will break, hence 9999
            if not self.fields_f["gene"] == self.fields_r["gene"]: 
                return 9999
            self.ampsize = self.fields_r["end"] - self.fields_f["start"]

        elif self.kit == "taqman":
            self.fields_p = get_fields(self.id_p)
            if not (self.fields_f["gene"] == self.fields_r["gene"] == self.fields_p["gene"]): 
                return 9999
            self.ampsize = self.fields_r["end"] - self.fields_f["start"]

        return self.ampsize

    def get_primerattrs(self):
        # get some generic information about
        lst_sensicals_f = [primer.dg for primer in self.primers_f if primer.is_sensical]
        lst_sensicals_r = [primer.dg for primer in self.primers_r if primer.is_sensical]
        self.meandg_f = round(sum(lst_sensicals_f) / len(lst_sensicals_f), 2)
        self.meandg_r = round(sum(lst_sensicals_r) / len(lst_sensicals_r), 2)

        lst_sensicals_f = [primer.tm for primer in self.primers_f if primer.is_sensical]
        lst_sensicals_r = [primer.tm_revc for primer in self.primers_r if primer.is_sensical]
        self.meantm_f = round(sum(lst_sensicals_f) / len(lst_sensicals_f), 1)
        self.meantm_r = round(sum(lst_sensicals_r) / len(lst_sensicals_r), 1)

        self.nambig_f = get_fields(self.id_f)["mm"]
        self.nambig_r = get_fields(self.id_r)["mm"]
        self.nambig_total = self.nambig_f + self.nambig_r

        self.n_sensicals = len(lst_sensicals_f) + len(lst_sensicals_r)
        self.n_primers = len(self.primers_f) + len(self.primers_r)

        if self.kit == "taqman":
            lst_sensicals_p = [primer.dg for primer in self.primers_p if primer.is_sensical]
            self.meandg_p = round(sum(lst_sensicals_p) / len(lst_sensicals_p), 2)

            lst_sensicals_p = [primer.tm for primer in self.primers_p if primer.is_sensical]
            self.meantm_p = round(sum(lst_sensicals_p) / len(lst_sensicals_p), 1)

            self.nambig_p = get_fields(self.id_p)["mm"]
            self.nambig_total += self.nambig_p
            self.n_sensicals += len(lst_sensicals_p)
            self.n_primers += len(self.primers_p)

    def get_alignment_info(self, internal_controls=[]):
        """ Given sets of forward, (probe), reverse primers, figure out how each variant of
        these primers align with each other, and if it doesn't work (i.e. bad gc/tm), return string
        if al conditions pass, return a list of all variant alignments """
        self.alignments = []
        lst_tmdiff_fr = []

        if self.kit == "endpoint":
            for primer_f in self.primers_f:
                for primer_r in self.primers_r:
                    flt_tmdiff = abs(primer_f.tm - primer_r.tm_revc)
                    flt_gcdiff = abs(primer_f.gc - primer_r.gc)

                    if primer_f.is_sensical and primer_r.is_sensical:
                        if flt_gcdiff > flt_maxgcdiff: 
                            return [], "fwd/rev GC difference too large"
                        if flt_tmdiff > flt_maxtmdiff: 
                            return [], "fwd/rev Tm difference too large"

                    aln_fr = get_alignments(primer_f.seq, revc(primer_r.seq, rev=False))
                    if aln_fr["dg"] < flt_mindg: 
                        return [], "fwd/rev dG too low"

                    # get internal control data
                    for rec_ic in internal_controls:
                        # all internal control sequences in file run 5>3
                        aln_fi = get_alignments(primer_f.seq, rec_ic.seq[::-1])
                        if aln_fi["dg"] < flt_mindg:
                            return [], "fwd interacts with IC"
                        aln_ri = get_alignments(primer_r.revc, rec_ic.seq[::-1])
                        if aln_ri["dg"] < flt_mindg:
                            return [], "rev interacts with IC"

                    aln_fr["pname_1"] = primer_f.id; aln_fr["pseq_1"] = primer_f.seq
                    aln_fr["pname_2"] = primer_r.id; aln_fr["pseq_2"] = primer_r.seq
                    aln_fr["tm_diff"] = flt_tmdiff
                    aln_fr["gc_diff"] = flt_gcdiff
                    dic_alnres = {"fr": aln_fr}
                    self.alignments.append(dic_alnres)
                    lst_tmdiff_fr.append(abs(primer_r.tm_revc - primer_f.tm))

            lst_dg = [aln["fr"]["dg"] for aln in self.alignments]
            self.meandg_fr = sum(lst_dg) / len(lst_dg)

        elif self.kit == "taqman":
            lst_dg = []; lst_tmdiff_p = []
            for primer_f in self.primers_f:
                for primer_r in self.primers_r:

                    # check forward/reverse combinations
                    flt_tmdiff_fr = abs(primer_r.tm - primer_f.tm)
                    flt_gcdiff_fr = abs(primer_r.gc - primer_f.gc)

                    if primer_f.is_sensical and primer_r.is_sensical: # "^" is python's xor operator / hz250
                        if flt_gcdiff_fr > flt_maxgcdiff:
                            return [], "fwd/rev GC difference too large"
                        if flt_tmdiff_fr > flt_maxtmdiff:
                            return [], "fwd/rev Tm difference too large"

                    # we also want to keep the primer names in the alignment result, so the function can return them
                    aln_fr = get_alignments(primer_f.seq, revc(primer_r.seq, rev=False))
                    aln_fr["pname_1"] = primer_f.id; aln_fr["pseq_1"] = primer_f.seq
                    aln_fr["pname_2"] = primer_r.id; aln_fr["pseq_2"] = primer_r.seq                    
                    aln_fr["tm_diff"] = flt_tmdiff_fr
                    aln_fr["gc_diff"] = flt_gcdiff_fr

                    if aln_fr["dg"] < flt_mindg:
                        return [], "fwd/rev dG too low"

                    # Forward/reverse combo OK - check probe
                    for primer_p in self.primers_p:
                        flt_probediff = primer_p.tm - (primer_f.tm + primer_r.tm) / 2
                        lst_tmdiff_p.append(flt_probediff)

                        # hz383 check internal controls - all three needs to be positive for set to be invalid
                        for rec_ic in internal_controls:
                            # all internal control sequences in file run 5>3
                            int_baddg = 0
                            aln_fi = get_alignments(primer_f.seq, rec_ic.seq[::-1])
                            if aln_fi["dg"] < flt_mindg: int_baddg += 1
                            aln_ri = get_alignments(primer_r.revc, rec_ic.seq[::-1])
                            if aln_ri["dg"] < flt_mindg: int_baddg += 1
                            aln_pi = get_alignments(primer_p.seq, rec_ic.seq[::-1])
                            if aln_pi["dg"] < flt_mindg: int_baddg += 1

                            if int_baddg >= 3:
                                if self.strand_p == "antisense":
                                    return [], "prb sense failed ~ antisense IC dG too low"
                                else:
                                    return [], "check antisense"
                            
                        # check for probe self-attributes if necessary
                        if self.strand_p == "antisense":
                            if flt_probediff < lst_tmprobes[0] or flt_probediff > lst_tmprobes[1]:
                                return [], "prb sense failed ~ antisense Tm not within range"

                        elif primer_p.is_sensical:
                            if flt_probediff < lst_tmprobes[0] or flt_probediff > lst_tmprobes[1]:
                                return [], "check antisense"

                        aln_fp = get_alignments(primer_f.seq, primer_p.seq[::-1])
                        aln_fp["pname_1"] = primer_f.id; aln_fp["pseq_1"] = primer_f.seq
                        aln_fp["pname_2"] = primer_p.id; aln_fp["pseq_2"] = primer_p.seq
                        aln_fp["gc_diff"] = "N/A"; aln_fp["tm_diff"] = "N/A"
                        if aln_fp["dg"] < flt_mindg:
                            if self.strand_p == "sense": return [], "check antisense"
                            else: return [], "prb sense failed ~ heterodimer antisense dG too low"
                            
                        aln_pr = get_alignments(primer_p.seq, revc(primer_r.seq, rev=False))
                        aln_pr["pname_1"] = primer_p.id; aln_pr["pseq_1"] = primer_p.seq
                        aln_pr["pname_2"] = primer_r.id; aln_pr["pseq_2"] = primer_r.seq
                        aln_pr["gc_diff"] = "N/A"; aln_pr["tm_diff"] = "N/A"
                        if aln_pr["dg"] < flt_mindg:
                            if self.strand_p == "sense": return [], "check antisense"
                            else: return [], "prb sense failed ~ heterodimer antisense dG too low"

                    # current variant combo passed requirements, move onto next
                    dic_alnres = {"fr": aln_fr, "fp": aln_fp, "pr": aln_pr}
                    lst_dg += [aln_fr["dg"], aln_fp["dg"], aln_pr["dg"]]
                    self.alignments.append(dic_alnres)
                    lst_tmdiff_fr.append(primer_r.tm_revc - primer_f.tm)

            # summarize combo stats
            #self.meantm_p = sum(lst_tmdiff_p) / len(lst_tmdiff_p)
            self.meandg_fpr = sum(lst_dg) / len(lst_dg)

        self.meantm_fr = sum(lst_tmdiff_fr) / len(lst_tmdiff_fr)

        return self.alignments, ""

# ----------------------------------------[ Generic Functions ]---------------------------------------

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

# 1. For Delta-G relative estimations
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

def write_alignment(wf, aln, p1_notation="sense", p2_notation="sense"):
    wf.write("Fwd ID/Rev ID: [{p1}/{p1}]\n".format(p1=aln["pname_1"], p2=aln["pname_2"]))
    wf.write("Fwd({n1})/Rev({n2}):\n[{s1}/{s2}]\n".format(n1=p1_notation, n2=p2_notation, s1=aln["pseq_1"], s2=aln["pseq_2"]))
    try: wf.write("Tm-diff: {0}\tdelta-G: {1}\n".format(round(aln["tm_diff"], 2), round(aln["dg"], 2)))
    except TypeError: wf.write("Tm-diff: {0}\tdelta-G: {1}\n".format(aln["tm_diff"], round(aln["dg"], 2)))
    wf.write("{0}\n\n".format(aln["printout"]))

# ======================================[ INITIALIZING PROGRAM ]======================================

# header block
print("[gs] Starting {0} for primer attribute analysis".format(str_program))
str_log = "\n--- STARTING [{0}] FOR PRIMER ATTRIBUTE ANALYSIS ---\n".format(str_program)
str_log += "\tPrimer input file:                            : {0}\n".format(pth_input.name)
str_log += "\tKit type:                                     : {0}\n".format(str_kit)
str_log += "\tLower and upper GC content limits             : {0}\n".format(lst_gclimits)
str_log += "\tLower and upper Tm limits                     : {0}\n".format(lst_tmlimits)
str_log += "\tLower and upper amplicon size limits          : {0}\n".format(lst_amplimits)
str_log += "\tCheck reference to reduce false negatives     : {0}\n".format(str_reference)
logger.info(str_log+"\n")

# ------------------------------------[ Getting input information ]-----------------------------------

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

# hz627 - getting internal control sequences
if str_kit == "endpoint":
    pth_ic = dir_ic / "IC_endpoint.fasta"
elif str_kit == "taqman":
    pth_ic = dir_ic / "IC_taqman.fasta"
with open(pth_ic) as handle:
    lst_icrecords = list(SeqIO.parse(handle, "fasta"))
print("[gs] Retrieved {0} internal control sequences from file: {1}".format(len(lst_icrecords), pth_ic.name))
logger.info("Retrieved {0} IC sequences from: {1}\n".format(len(lst_icrecords), pth_ic))


# Getting primer information into dictionary
print("[gs] Indexing primers from file: {0}".format(pth_input.name))
# Read primer input file into dictionary
with open(pth_input, "rt") as rf_input:
    lst_lines = rf_input.readlines()

current_primers = {}
pname_prev = ""
for i in range(0, len(lst_lines), 2):
    pname_exp = lst_lines[i].strip()[1:]
    pname_com = shorten_name(pname_exp)
    primer = Primer(pname_exp, lst_lines[i+1].strip())

    if pname_com == pname_prev: current_primers[pname_prev].append(primer)
    else: current_primers[pname_com] = [primer]
    pname_prev = pname_com

# remove unused variables
del lst_lines

# ========================================[ FILTRATION STEPS ]========================================
int_initial = len(current_primers)
logger.info("Starting number of primers (combined)          : {0}\n".format(int_initial))
print("\n---> {0} Starting primers prior to filtration".format(int_initial))

# ---------------------------[ Getting reference from file (if applicable) ]--------------------------
if str_reference != "none":
    gene_seqs = {}
    pth_reference = dir_project / str_reference
    print("[gs] Reading reference genome into program: {0}".format(pth_reference.name))
    with open(pth_reference) as handle:
        lst_reference = list(SeqIO.parse(handle, "fasta"))

    for record in lst_reference:
        gene_seqs[record.id] = record.seq.strip().replace("-","").upper()
    del lst_reference

    for pname_com, lst_primers in current_primers.items():
        for primer in lst_primers:
            primer.check_reference()

# -------------------------------------[ Calculating GC content ]-------------------------------------

if 'gc' in lst_runtype:
    print("[gs] Filtering primers by gc content: {0}".format(lst_gclimits))
    logger.info("Filtering primers by GC content...")

    lst_removals = []
    for pname_com, lst_primers in current_primers.items():
        for primer in lst_primers:
            flt_gc = primer.get_gc()
            # if the GC content is good, we nove onto the next primer variant
            if flt_gc >= lst_gclimits[0] and flt_gc <= lst_gclimits[1]: continue
            # break means we're done with this primer and all its variants - nothing more to do
            if primer.is_sensical: lst_removals.append(pname_com); break

    for id_com in lst_removals:
        current_primers.pop(id_com, None)

    print("---> {0} primers removed for GC content".format(len(lst_removals)))
    logger.info("Primers (combined) removed by GC content: {0}".format(len(current_primers)))


# ---------------------------------[ Estimating Melting Temperature ]---------------------------------

if 'tm' in lst_runtype:
    print("[gs] Filtering primers by melting temperature: {0}".format(lst_tmlimits))
    logger.info("Filtering primers by Tm values...")

    lst_removals = []
    for pname_com, lst_primers in current_primers.items():
        for primer in lst_primers:
            flt_tm, flt_tm_revc = primer.get_tm()
            # if the Tm is good (either the sense or antisense), we move onto the next primer variant
            if max(flt_tm, flt_tm_revc) >= lst_tmlimits[0] and min(flt_tm, flt_tm_revc) <= lst_gclimits[1]: 
                continue
            if primer.is_sensical:
                lst_removals.append(pname_com); break
    
    for id_com in lst_removals:
        current_primers.pop(id_com, None)

    print("---> {0} primers removed for Tm values".format(len(lst_removals)))
    logger.info("Primers (combined) removed by Tm values: {0}".format(len(current_primers)))


# --------------------------------[ Estimating Delta-G for homodimers ]-------------------------------

if 'dg' in lst_runtype:
    print("[gs] Estimating binding strength for remaining primers")
    logger.info("Filtering primers by homodimer stability...")

    lst_removals = []
    for pname_com, lst_primers in current_primers.items():
        for primer in lst_primers:
            flt_dg = primer.get_dg()
            if flt_dg >= flt_mindg: continue
            if primer.is_sensical: 
                lst_removals.append(pname_com); break

    for id_com in lst_removals:
        current_primers.pop(id_com, None)

    print("---> {0} primers removed for homodimers.".format(len(lst_removals)))

print("---> {0} of {1} primers remaining after preliminary filters.\n".format(len(current_primers), int_initial))
logger.info("Primers remaining before generating combinations: {0}\n".format(len(current_primers)))


# ----------------------------------[ Write current primers to file ]---------------------------------

print("[gs] Writing sequences to file: {0}".format(pth_out_primers.name))
with open(pth_out_primers, "wt") as wf_primers:
    for pname_com, lst_primers in current_primers.items():
        for primer in lst_primers:
            wf_primers.write(">{0}\n".format(primer.id))
            wf_primers.write("{0}\n".format(primer.seq))

if tf_primersonly:
    print("[gs] Program complete.")
    quit()


# =====================================[ GENERATING PRIMER SETS ]=====================================

logger.info("Filtering primers by amplicon size and heterodimer stability")
# 2D sets for primer combinations
lst_removals = [] # primer combined names to remove from pool
dic_removals = {} # reason to remove primer from pool (for diagnostic purposes)
dic_combos = {} # [(primer1.id, primer2.id)] = Pcombo()

if str_kit == "endpoint":
    # ----------------------------------[ Endpoint Primer Sets ]--------------------------------------

    print("[gs] Gathering amplicon combinations for EP kits")
    # First, we want to get the possible primers that will give us an appropriate amplicon

    for p1, (pname_com1, lst_primers1) in enumerate(current_primers.items()):
        for p2, (pname_com2, lst_primers2) in enumerate(current_primers.items()):
            if p1 >= p2: continue

            combo = Pcombo(pname_com1, pname_com2, lst_primers1, lst_primers2)
            int_ampsize = combo.get_ampsize()
            if int_ampsize > lst_amplimits[1]: break
            elif int_ampsize < lst_amplimits[0]: continue
            dic_combos[(pname_com1, pname_com2)] = combo

    int_initialcombos = len(dic_combos)
    print("---> {0} possible number of F/R combinations identified.\n".format(int_initialcombos))
    logger.info("Possible fwd/rev combinations identified: {0}".format(int_initialcombos))

    # pick out which primer sets to remove
    del current_primers
    print("[gs] Filtering endpoint primer sets based on biochemical properties...")

    with open(pth_out_aln, "wt") as wf_aln, open(pth_out_psets, "wt") as wf_psets, open(pth_out_fails, "wt") as wf_pfails:
        print("[gs] Writing primer alignments to file: {0}".format(pth_out_aln.name))
        print("[gs] Writing primer combinations to file: {0}".format(pth_out_psets.name))

        # make sure to not use special characters because these headers will be used as class attributes later
        wf_psets.write("id_f,tm_f,dg_f,id_r,tm_r,dg_r,n_primers,amp_size,tm_diff_fr,dg_het\n")
        wf_pfails.write("id_f,tm_f,dg_f,id_r,tm_r,dg_r,n_primers,amp_size,reason\n")

        for i, (key, combo) in enumerate(dic_combos.items()):
            lst_alignments, msg = dic_combos[key].get_alignment_info(internal_controls=lst_icrecords)
            combo.get_primerattrs()

            if len(lst_alignments) < 1:
                if msg in dic_removals: dic_removals[msg] += 1
                else: dic_removals[msg] = 1
                lst_removals.append(key)
                wf_pfails.write("{id},{tm},{dg},".format(id=combo.id_f, tm=combo.meantm_f, dg=combo.meandg_f))
                wf_pfails.write("{id},{tm},{dg},".format(id=combo.id_r, tm=combo.meantm_r, dg=combo.meandg_r))
                wf_pfails.write("{n},{amp},{msg}\n".format(n=combo.n_primers, amp=combo.ampsize, msg=msg))
                
            else:
                # write combo to file
                wf_psets.write("{id},{tm},{dg},".format(id=combo.id_f, tm=combo.meantm_f, dg=combo.meandg_f))
                wf_psets.write("{id},{tm},{dg},".format(id=combo.id_r, tm=combo.meantm_r, dg=combo.meandg_r))
                wf_psets.write("{n},{amp},{tm},{dg}\n".format(n=combo.n_primers, amp=combo.ampsize, tm=combo.meantm_fr, dg=combo.meandg_fr))
                wf_aln.write("PCombo #{0}\n".format(i))
                for aln_set in lst_alignments:
                    write_alignment(wf_aln, aln_set["fr"], p1_notation="sense", p2_notation="sense")

    for key in lst_removals:
        dic_combos.pop(key, None)

elif str_kit == "taqman":
    # -----------------------------------[ taqman Primer Sets ]--------------------------------------

    print("[gs] Gathering amplicon combinations for TM kits")
    for p1, (pname_com1, lst_primers1) in enumerate(current_primers.items()):
        for p2, (pname_com2, lst_primers2) in enumerate(current_primers.items()):
            if p2 <= p1: continue
            # temporarily create a Pcombo class to get ampsize, so we can skip to next p1 if too large
            _combo = Pcombo(pname_com1, pname_com2, lst_primers1, lst_primers2)
            int_ampsize = _combo.get_ampsize()
            if int_ampsize > lst_amplimits[1]: break
            elif int_ampsize < lst_amplimits[0]: continue
            
            #hz838 TODO - check if probe starts with a 'G' or has four consecutive 'G's in it - we don't want that
            del _combo, int_ampsize


            
            for p3, (pname_com3, lst_primers3) in enumerate(current_primers.items()):
                if p3 <= p2: continue
                # check that the end of the probe starts before the start of the reverse
                elif get_fields(pname_com3)["start"] - get_fields(pname_com2)["end"] < 0: continue

                # p1:forward, p2:probe, p3:reverse
                combo = Pcombo(pname_com1, pname_com3, lst_primers1, lst_primers3, kit=str_kit)
                combo.add_probe(pname_com2, lst_primers2)
                int_ampsize = combo.get_ampsize()
                if int_ampsize > lst_amplimits[1]: break
                elif int_ampsize < lst_amplimits[0]: continue
                dic_combos[(pname_com1, pname_com2, pname_com3)] = combo

    # pick out primers to remove
    del current_primers
    int_initialcombos = len(dic_combos)
    print("---> {0} possible number of F/P/R combinations identified.\n".format(int_initialcombos))
    logger.info("Possible fwd/prb/rev combinations identified: {0} ".format(int_initialcombos))

    with open(pth_out_aln, "wt") as wf_aln, open(pth_out_psets, "wt") as wf_psets, open(pth_out_fails, "wt") as wf_pfails:
        print("[gs] Writing primer alignments to file: {0}".format(pth_out_aln.name))
        print("[gs] Writing primer combinations to file: {0}".format(pth_out_psets.name))

        wf_psets.write("id_f,tm_f,dg_f,id_p,strand_p,tm_p,dg_p,id_r,tm_r,dg_r,n_primers,amp_size,tm_diff_fr,dg_het\n")
        wf_pfails.write("id_f,tm_f,dg_f,id_p,tm_p,dg_p,id_r,tm_r,dg_r,n_primers,amp_size,reason\n")

        # pick out which primer sets to remove
        print("[gs] Filtering taqman primer sets based on biochemical properties...")
        for i, (key, combo) in enumerate(dic_combos.items()):
            # we need to compare three alignments
            lst_alignments, msg = dic_combos[key].get_alignment_info(internal_controls=lst_icrecords)
            combo.get_primerattrs()

            # if the probe failed, we can try the antisense version to see if we can salvage the combo
            if msg == "check antisense":
                # get the complement of the 5>3s for all variants of the probe
                primers_p = [Primer(primer.id, revc(primer.seq, rev=False)) for primer in combo.primers_p]

                # check that the new antisense probe meets GC, Tm, dG requirements
                int_pchecked = 0
                for primer in primers_p:
                    flt_gc = primer.get_gc()
                    flt_tm, flt_tm_revc = primer.get_tm()
                    flt_dg = get_alignments(primer.seq, primer.seq[::-1])["dg"]

                    if not primer.is_sensical:
                        int_pchecked += 1
                        continue # we're not checking self-attributes for non-sensical probes

                    if flt_gc < lst_gclimits[0] or flt_gc > lst_gclimits[1]: break
                    elif max(flt_tm, flt_tm_revc) < lst_tmlimits[0] or min(flt_tm, flt_tm_revc) > lst_tmlimits[1]: break                    
                    elif flt_dg < flt_mindg: break
                    int_pchecked += 1
                
                if int_pchecked == len(primers_p):
                    dic_combos[key].add_probe(combo.id_p, primers_p, strand_p="antisense")
                    lst_alignments, msg = dic_combos[key].get_alignment_info(internal_controls=lst_icrecords)
                else:
                    lst_alignments = []
                    msg = "prb sense failed ~ antisense did not meet self Tm/dG requirements"

            if len(lst_alignments) < 1:
                if msg in dic_removals: dic_removals[msg] += 1
                else: dic_removals[msg] = 1
                lst_removals.append(key)

                wf_pfails.write("{id},{tm:.1f},{dg:.2f},".format(id=combo.id_f, tm=combo.meantm_f, dg=combo.meandg_f))
                wf_pfails.write("{id},{tm:.1f},{dg:.2f},".format(id=combo.id_p, tm=combo.meantm_p, dg=combo.meandg_p))
                wf_pfails.write("{id},{tm:.1f},{dg:.2f},".format(id=combo.id_r, tm=combo.meantm_r, dg=combo.meandg_r))
                wf_pfails.write("{n},{amp},{msg}\n".format(n=combo.n_primers, amp=combo.ampsize, msg=msg))

            else:
                wf_psets.write("{id},{tm:.1f},{dg:.2f},".format(id=combo.id_f, tm=combo.meantm_f, dg=combo.meandg_f))
                wf_psets.write("{id},{s},{tm:.1f},{dg:.2f},".format(id=combo.id_p, s=combo.strand_p, tm=combo.meantm_p, dg=combo.meandg_p))
                wf_psets.write("{id},{tm:.1f},{dg:.2f},".format(id=combo.id_r, tm=combo.meantm_r, dg=combo.meandg_r))
                wf_psets.write("{n},{amp},".format(n=combo.n_primers, amp=combo.ampsize))
                wf_psets.write("{tm:.1f},{dg:.2f}\n".format(tm=combo.meantm_fr, dg=combo.meandg_fpr))

                # hz916 not writing the alignment to save space
                # wf_aln.write("PCombo #{0}\n:".format(i))
                # for aln_set in lst_alignments:
                #     write_alignment(wf_aln, aln_set["fr"])
                #     write_alignment(wf_aln, aln_set["fp"], p2_notation=combo.strand_p)
                #     write_alignment(wf_aln, aln_set["pr"], p1_notation=combo.strand_p)
        
    # we don't really have to do this.. 
    # just thought it'd be easier in case we have to do some work with combos afterwards
    for key in lst_removals:
        dic_combos.pop(key, None)

# =============================================[ FOOTER ]=============================================

print("---> {0} combinations removed for heterodimer formation.".format(len(lst_removals)))
print("---> {0} of {1} combinations remaining after filtration.".format(len(dic_combos), int_initialcombos))

str_log = "Primers removed:\n"
for msg, count in dic_removals.items():
    str_log += "\t{0}: {1}\n".format(msg, count)
str_log += "Combos written: {0}\n".format(len(dic_combos))
str_log += "\nProgram complete.\n"+"-"*80+"\n"
logger.info(str_log)

print("[gs] Program complete.")
