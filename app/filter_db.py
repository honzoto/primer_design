#!/usr/bin/python3

import sys, os, re, glob
import textwrap
from pathlib import Path
import pandas as pd
import logging

from Bio import SeqIO

"""
PURPOSE

The objective is to find primers specific to one group (positive) and does not align to another (negative)
We've generated a list of primers (fasta format) based on the positive dataset, and have aligned the list
to a BLAST databse creative from the negative dataset. We now need to generate a filtered primer list,
filtering out the primers that aligned to the negatives.

- FASTA files have to be called *.fasta, NOT *.fa or *.fna

making blastdb command:
makeblastdb -in High_risk_strains.fasta -title "High_risk_HPV" -dbtype 'nucl' -logfile log_hpv.txt

sample BLAST search through CL:
blastn -db /mnt/tank/bench/projects/Primer_design/Highrisk_database/High_risk_strains.fasta 
-query primerlist_group1_alignments_3mm_exp.fasta 
-out primerlist_group1_alignments_3mm_exp_blastoutput.txt 

UPDATE NOTES
v1  - we can now perform the search in this program too in addition to filtering
    - the filtering algorithm now uses pandas to read the blast result as a dataframe, instead of parsing
        the text file line-by-line
    - added i/o params as an alternative way of passing primer file into program (for pipeline compatibility)
    - working on option to give user the ability to select multiple negative databases
v2  - using logging library in place of manual logging
    - update to remove_primers function, internal control removal structure
v3  - adding positive database search at the end of negative database filtering
    - new method of parsing blast search result (line-by-line) - doing so manually to get primer counts
    - primers removed metric does not actually measure primers removed, (bc some primers are not in dic_primers_exp)
        - we fixed that by outputting int_removed from the removal function
    - fixed up some bugs with the positive database search
    - moving database search results into its own folder inside the project
v4  - script will now perform positive db search before negative
    - a temporary FASTA file will be created, and reduced each time primers are eliminated
    - removed internal control filtering (moved to generate_psets.py)
    - hardcoded pct coverage (set to 80), evalue to 1000
---
- we haven't used the global alignment algorithm with BLAST search that we're using for covid_primers
    - we may not even have to implement this unless we're not getting enough primers at the end
- we don't have the same primer class setup we have with generate_psets.py, which makes handling variants more organized
> we're using the entire FASTA file as a query for blastn every time, however, after removing some sequences,
    we could just write into a temporary fasta file and use that for subsequent searches - just a thought

"""

# =======================================[ DECLARING VARIABLES ]======================================

dir_program = Path("/mnt/tank/bench/scripts/primer_design")
str_version = "1.4.2"
flt_minpct = 0.80 # has to hit at least this pct of sequences in positive DB
str_output = "auto"
str_searchparams = " -num_threads 20 -outfmt 7 -task megablast -word_size 16 -qcov_hsp_perc 100 -evalue 1000"
lst_negdb = []; lst_posdb = []

# Get a list of all available database folders , putting them into string format
os.chdir(dir_program)
dbfolders = sorted([Path(p).name for p in glob.glob("database/*")])
str_dbfolders = "\n    ".join(textwrap.wrap(", ".join(dbfolders)))

str_helpmenu = \
"""
----------------------------------[ HELP MENU ]---------------------------------

    -h/--help       : shows this menu
    -p/--project    : <dir> name of the project folder
    -i/--input      : <path> name of the FASTA file from running the previous 
                        script: primer_figurator.py
    -o/--output     : <path> name of output FASTA file
    -x/--threshold  : <float> [default={x}] primers must align to at least 
                        this percent of sequences in positive databases
    -pd/--pos_db    : <dir> name of positive db folder(s) in 'database' 
                        directory - only keep primers that align
    -nd/--neg_db    : <dir> name of negative db folder(s) in 'database' 
                        directory - filter out primers that align

    Available databases:
    {db}

--------------------------------------------------------------------------------
""".format(x=flt_minpct, db=str_dbfolders)

# -----------------------------------[ Getting arguments from user ]----------------------------------
if len(sys.argv) > 1:
    for i, arg in enumerate(sys.argv):
        if arg == ("-p" or "--project"):
            str_project = sys.argv[i+1]
        elif arg == ("-i" or "--input"):
            str_primerfile = sys.argv[i+1]
        elif arg == ("-o" or "--output"):
            str_output = sys.argv[i+1]
        elif arg == ("-x" or "--threshold"):
            flt_minpct = float(sys.argv[i+1])
        elif arg == ("-pd" or "--pos_db"):
            lst_posdb = sys.argv[i+1].strip().split(",")
        elif arg == ("-nd" or "--neg_db"):
            lst_negdb = sys.argv[i+1].strip().split(",")
        elif arg == ("-h" or "--help"):
            print(str_helpmenu); quit()
        elif arg[0] == "-":
            print("[WARNING] cannot recognize argument: {0}".format(arg))
else:
    print("[ERROR] no arguments specified.")
    print(str_helpmenu); quit()

# -----------------------------[ Setting directories and user variables ]-----------------------------

dir_project = Path("/mnt/tank/bench/projects") / str_project
pth_input = dir_project / str_primerfile
pth_log = dir_project / "PrimerDesign.log"

if str_output == "auto": pth_outfile = dir_project / (str_primerfile.split(".")[0]+"_filtered.fasta")
else: pth_outfile = dir_project / str_output

# Setting up logger
print("[fd] Setting up logger: {0}".format(pth_log))
#open(pth_log, "wt").close() # reset the log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
log_formatter = logging.Formatter("%(asctime)s-%(levelname)s: %(message)s", datefmt="[%Y-%m-%d %H:%M:%S]")
file_handler = logging.FileHandler(pth_log)
file_handler.setFormatter(log_formatter)
logger.addHandler(file_handler)


# ======================================[ CLASSES AND FUNCTIONS ]=====================================

def hamming_distance(chain1, chain2):
    return sum(c1 != c2 for c1, c2 in zip(chain1, chain2))

def shorten_name(expanded_name):
    try: # see if there's an easier way to do this regex search
        m = re.search('^>?([^\|]*)\|([^\|]*)\|(\d+)\|(\d+)(\~\d*)?', expanded_name)
        return "{0}|{1}|{2}|{3}".format(m.group(1), m.group(2), m.group(3), m.group(4))
    except:
        logger.warning("Cannot shorten primer name: {0}".format(expanded_name))
        return expanded_name

def remove_primers(lst_names_com, dic_primers_exp):
    """ Accepts a list of shortened primer names (e.g. E7|108-131|24|1)
    and removes entries that matches the pattern (e.g. E7|108-131|24|1~1, E7|108-131|24|1~2) 
    - i don't think we have to return the same dictionary, because the changes are kept
    """
    lst_removals = [name for name in dic_primers_exp.keys() if shorten_name(name) in lst_names_com]
    for str_name_exp in lst_removals:
        dic_primers_exp.pop(str_name_exp, None)

    set_removals = set([shorten_name(name) for name in lst_removals])
    return len(set_removals)

def check_match(blastres, int_maxmm=1):
    # Read the BLAST result as a dataframe and filter down based on mismatches, gaps
    lst_header_ = "query acc.ver,subject acc.ver,% identity,alignment length,mismatches,\
        gap opens,q. start,q. end,s. start,s. end,evalue,bit score".split(",")
    lst_header = [v.strip() for v in lst_header_]

    df_blastres = pd.read_table(blastres, comment="#", names=lst_header)
    df_subset = df_blastres[df_blastres['gap opens'] == 0]
    df_subset = df_subset[df_subset['mismatches'] <= int_maxmm]
    
    # we don't care about one primer hitting two negative-db accessions, remove it regardless
    lst_subset = list(set(df_subset['query acc.ver']))
    return lst_subset

def check_alignment(query, reference, allowed_mm=0):
    """ check if you can find the query in a larger reference sequence,
    the allowed mm is the maximum hamming distance"""

    if allowed_mm < 1: # if not in reference, -1 will be returned
        return reference.find(query)

    # if we allow mismatches, then we have to check each nucleotide in the sequence
    for r in range(len(reference) - len(query)):
        subref = reference[r:r+len(query)]
        if hamming_distance(query, subref) < allowed_mm:
            return r
    return -1

# ==========================================[ PROGRAM START ]=========================================

# header block
print("[fd] Starting blast_filterer.py")
str_log = "\n--- STARTING blast_filterer.py ---\n"
str_log += "\tPrimer input file     : {0}\n".format(pth_input)
str_log += "\tNeg DB name(s)        : {0}\n".format(lst_negdb)
str_log += "\tPos DB name(s)        : {0}\n".format(lst_posdb)
str_log += "\tPos min. db coverage  : {0}\n".format(flt_minpct)
logger.info(str_log+"\n")

# create a temporary copy of the input file to use in workflow

print("[fd] Indexing primers as SeqIO records")
dic_primers_exp = SeqIO.to_dict(SeqIO.parse(pth_input, 'fasta')) # if the above line doesn't work
int_initial_exp = len(dic_primers_exp)
logger.info("Starting number of oligos (expanded): {0}".format(int_initial_exp))

# Check number of non-ambiguous primers in pool
set_primers_com = set([shorten_name(k) for k in dic_primers_exp.keys()])
int_initial_com = len(set_primers_com)
logger.info("Starting number of primers (combined): {0}\n".format(len(set_primers_com)))
print("---> {0} starting number of primers".format(len(set_primers_com)))

# first, create a folder in the projects to keep all our search results
try: os.mkdir(dir_project / "search_results")
except: print("[fd] search_results folder already exists.")

# =================================[ FILTERING BY POSITIVE DATABASES ]================================

str_logremoval = "Primers removed by positive databases:\n"
for str_posdb in lst_posdb:
    str_blastres = "blastres_{0}.txt".format(str_posdb)
    pth_blastres = dir_project / "search_results" / str_blastres

    # ----------------------------------[ Perform BLAST search ]--------------------------------------
    print("\n[fd] Searching positive database: {0}".format(str_posdb))
    str_dbre = "/mnt/tank/bench/scripts/primer_design/database/{0}/*.fasta".format(str_posdb)
    pth_fastadb = glob.glob(str_dbre)[0] # hopefully there's only one .fasta file in there
    
    try:
        #with open(pth_makedblog, "rt") as rf_dblog:
        str_dbre = "/mnt/tank/bench/scripts/primer_design/database/{0}/log*.txt".format(str_posdb)
        pth_makedblog = glob.glob(str_dbre)[0]
        with open(pth_makedblog, "rt") as rf_dblog:
            for line in rf_dblog:
                if line.find("Adding sequences from FASTA;") >= 0:
                    int_numseqs = int(line.split()[5]); break
    except:
        # if we didn't make a log during DB creation, we have to manually count the number of lines in the FASTA
        print("[WARNING] no log file found for database creation. Counting FASTA lines from file instead...")
        logger.warning("No log file found for positive DB creation... counting manually.".format(str_posdb))
        int_numseqs = len([1 for line in open(pth_fastadb, "r") if line.startswith(">")])
    
    logger.info("db path: {0}\n\tsequences found: {1}".format(pth_fastadb, int_numseqs))
    print("{0} total sequences found in db: {1}".format(int_numseqs, str_posdb))

    # search BLAST database with primer output from previous script
    str_search = "blastn -db {db} -query {q} -out {o} ".format(db=pth_fastadb, q=pth_input, o=pth_blastres) 
    print("[cmd] " + str_search + str_searchparams)
    logger.info("[cmd] {0}{1}".format(str_search, str_searchparams))
    os.system(str_search + str_searchparams)

    # --------------------------[ Filtering primers with search output ]------------------------------

    dic_primerhits = {}
    print("[fd] reading blast output from file:", pth_blastres)

    with open(pth_blastres, "rt") as rf_blastres:
        for line in rf_blastres:
            if line.startswith("#"): continue
            lst_line = line.strip().split("\t")
            str_pname_com = shorten_name(lst_line[0])
            if str_pname_com in dic_primerhits:
                dic_primerhits[str_pname_com].add(lst_line[1])
            else: # watch out when creating 1-element sets, because set('abc') >>> {'a', 'b', 'c'}
                dic_primerhits[str_pname_com] = set([lst_line[1]])
                
    set_removals_com = set()
    for str_pname_com in dic_primerhits:
        if len(dic_primerhits[str_pname_com]) / int_numseqs < flt_minpct:
            #print("hz308", len(dic_primerhits[str_pname_com]) / int_numseqs, str_pname_com)
            set_removals_com.add(str_pname_com)

    int_removed = remove_primers(set_removals_com, dic_primers_exp)
    str_logremoval += "\t[{0}]: {1}\n".format(str_posdb, int_removed)
    print("---> {0} primers removed from pool by pos-db: {1}".format(int_removed, str_posdb))

logger.info(str_logremoval)

# =================================[ FILTERING BY NEGATIVE DATABASES ]================================

# we do negative before positive because we want to cut down on the primer pool first before doing more work
str_logremoval = "Primers removed by negative databases:\n"
for str_negdb in lst_negdb:
    str_blastres = "blastres_{0}.txt".format(str_negdb)
    pth_blastres = dir_project / "search_results" / str_blastres

    # ----------------------------------[ Perform BLAST search ]--------------------------------------
    print("\n[fd] Searching negative database: {0}".format(str_negdb))
    str_dbre = "/mnt/tank/bench/scripts/primer_design/database/{0}/*.fasta".format(str_negdb)
    pth_fastadb = glob.glob(str_dbre)[0] # hopefully there's only one .fasta file in there    
    logger.info("db path: {0}\n".format(pth_fastadb))

    # search BLAST database with primer output from previous script
    str_search = "blastn -db {db} -query {q} -out {o} ".format(db=pth_fastadb, q=pth_input, o=pth_blastres) 
    print("[cmd] " + str_search + str_searchparams)
    logger.info("[cmd] {0}{1}".format(str_search, str_searchparams))
    os.system(str_search + str_searchparams)

    # --------------------------[ Filtering primers with search output ]------------------------------
    print("[fd] reading blast output from file:", pth_blastres)

    set_accessions_com = set([shorten_name(acc) for acc in check_match(pth_blastres)])
    int_removed = remove_primers(set_accessions_com, dic_primers_exp)
    str_logremoval += "\t[{0}]: {1}\n".format(str_negdb, int_removed)
    print("---> {0} primers removed from pool by neg-db: {1}".format(int_removed, str_negdb))

set_primers_com = set([shorten_name(k) for k in dic_primers_exp.keys()])
print("---> {0} primers remaining after negative-db filters".format(len(set_primers_com)))
str_logremoval += "\tPrimers remaining: {0} of {1}".format(len(set_primers_com), int_initial_com)
logger.info(str_logremoval)

# ================================[ WRITING REMAINING PRIMERS TO FILE ]===============================

set_primers_com = set([shorten_name(k) for k in dic_primers_exp.keys()])
int_final = len(set_primers_com)

print("\n[fd] Finished filtration steps...")
print("---> {0} primers remaining. Writing to file: {1}".format(int_final, pth_outfile))
with open(pth_outfile, 'w') as handle:
    SeqIO.write(dic_primers_exp.values(), handle, 'fasta')

logger.info("Program complete.\n"+"-"*80+"\n")
print("---> {f} of {i} primers passed filtration".format(f=int_final, i=int_initial_com))
print("[fd] program complete.")