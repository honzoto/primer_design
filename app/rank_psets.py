#!/usr/bin/python3

import sys, re
import pandas as pd
from pathlib import Path
import logging
from Bio import SeqIO

"""
PURPOSE:
After obtaining the psets output file from generate_psets.py, we want to know which primer sets are the best.
To do so, we can assign weights to individual attributes of each primer set and provide a score.
This score can be used to rank the sets from best to worst.

UPDATE NOTES:
v1  - program version for all scripts in the pipeline changed to 3-field notation
    - created processes to calculate score and sort set quality
    - created processes to retrieve original ambig primer sequences from ambig file
    - hz123 - auto-detect kit type based on number of columns in CSV file

How to sort dictionaries of objects by attribute value?
for student in (sorted(student_Dict.values(), key=operator.attrgetter('age'))):
    print(student.name)

"""
str_version = "1.1.2"
str_program = Path(__file__).name

str_kit = "auto"
str_primerfile = "none"

help_menu = """
----------------------------------[ HELP MENU ]---------------------------------
    
    USAGE: rank_psets.py -p <project> -i <input> -m <primers>

    < required arguments >
    -p/--project    : <directory> name of the project directory
    -i/--input      : <file> name of primer sets CSV from generate_psets.py

    < optional arguments >
    -r/--primers    : [default={m}] name of ambiguous FASTA file from 
                    get_primers.py (use non-ambiguous FASTA for if x=0 used)
    -k/--kit        : [default={k}] kit type used to generate primer sets
                        (options: auto/endpoint/taqman)
    -h/--help       : shows this menu

--------------------------------------------------------------------------------
""".format(m=str_primerfile, k=str_kit)


# ======================================[ VARIABLE DECLARATIONS ]=====================================

if len(sys.argv) > 1:
    lst_args = sys.argv
    for i in range(len(lst_args)):
        if lst_args[i] == ("-p" or "--project"):
            str_project = lst_args[i+1]
        elif lst_args[i] == ("-i" or "--input"):
            str_input = lst_args[i+1]
        elif lst_args[i] == ("-r" or "--primers"):
            str_primerfile = lst_args[i+1]
        elif lst_args[i] == ("-k" or "--kit"):
            str_kit = lst_args[i+1].lower()
        elif lst_args[i] == ("-h" or "--help"):
            print(help_menu); quit()

        elif lst_args[i][0] == "-":
            print("[WARNING] unrecognized flag: {0}".format(lst_args[i]))
else:
    print("Use -h or --help to see program options...\n")
    print(help_menu); quit()

dir_project = Path("/mnt/tank/bench/projects") / str_project
pth_input = dir_project / str_input
pth_output = dir_project / str_input.replace(".csv", "_sorted.csv")
pth_primers = dir_project / str_primerfile

# setting up logger
pth_log = dir_project / "PrimerDesign.log"
print("[de] Setting up logger: {0}".format(pth_log))
#open(pth_log, "wt").close() # reset the log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
log_formatter = logging.Formatter("%(asctime)s-%(levelname)s: %(message)s", datefmt="[%Y-%m-%d %H:%M:%S]")
file_handler = logging.FileHandler(pth_log)
file_handler.setFormatter(log_formatter)
logger.addHandler(file_handler)


# ======================================[ CLASSES AND FUNCTIONS ]=====================================

class Target:
    def __init__(self, labels: list, target=0.0, weight=0.0):
        self.labels = labels
        self.target = target
        self.weight = weight

def get_score(obs, tgt, wt=-1):
    if wt == -1: wt = tgt.weight
    try: return abs(round((obs - tgt.target) / tgt.target * wt, 2))
    except ZeroDivisionError: return abs(round((obs - tgt.target) / (tgt.target + 1) * wt, 2))

def shorten_name(expanded_name):
    try: 
        m = re.search('^>?([^\|]*)\|([^\|]*)\|(\d+)\|(\d+)(\~\d*)?', expanded_name)
        return "{0}|{1}|{2}|{3}".format(m.group(1), m.group(2), m.group(3), m.group(4))
    except:
        logger.warning("Cannot shorten primer name: {0}".format(expanded_name))
        return expanded_name

# =========================================[ PROCESSING DATA ]========================================

global tgt_nprimers, tgt_homodimer, tgt_heterodimer, tgt_tmdiff
tgt_nprimers = Target(["n_primers"], 0.0, 0.40)
tgt_tmdiff = Target(["tm_diff_fr"], 0.0, 0.25)
tgt_homodimer = Target(["dg_f", "dg_r", "dg_p"], 0.0, 0.20)
tgt_heterodimer = Target(["dg_fr", "dg_het"], 0.0, 0.15)

print("[rp] Using input file for target/adjustment determination: {0}".format(pth_input))
logger.info("Starting honzo's {0} v{1} for sorting primer sets.".format(str_program, str_version))
df_combos = pd.read_csv(pth_input)
if str_kit == "auto": # hz123 determine which kit type was used by getting number of columns
    if len(df_combos.columns) <= 10: str_kit = "endpoint"
    else: str_kit = "taqman"
logger.info("Kit type selected: {0}".format(str_kit))

# determining target values for attributes of interest
tgt_nprimers.target = min(df_combos["n_primers"])
tgt_homodimer.target = max(list(df_combos["dg_f"]) + list(df_combos["dg_r"]))
tgt_heterodimer.target = max(df_combos["dg_het"])
try: tgt_tmdiff.target = min(df_combos["tm_diff_fr"])
except: tgt_tmdiff.target = min(df_combos["tm_diff"])

# Writing target values to log file
str_log = "[Targets, weight]\n"
str_log += "\tNumber of primers:    {0:.2f}\twt={1:.2f}\n".format(tgt_nprimers.target, tgt_nprimers.weight)
str_log += "\tTm difference (F/R):  {0:.2f}\twt={1:.2f}\n".format(tgt_tmdiff.target, tgt_tmdiff.weight)
str_log += "\tdG homodimer:         {0:.2f}\twt={1:.2f}\n".format(tgt_homodimer.target, tgt_homodimer.weight)
str_log += "\tdG heterodimer:       {0:.2f}\twt={1:.2f}\n".format(tgt_heterodimer.target, tgt_heterodimer.weight)
logger.info(str_log+"\n")

# ---------------------------------[ Getting scores for primer sets ]---------------------------------

print("[rp] Getting scores for primer sets")
logger.info("Getting scores of primer sets by kit: {0}".format(str_kit))
ser_nprimers = get_score(df_combos["n_primers"], tgt_nprimers)
ser_tmdiff = get_score(df_combos["tm_diff_fr"], tgt_tmdiff)
if str_kit == "endpoint": 
    ser_homodimer_f = get_score(df_combos["dg_f"], tgt_homodimer, wt=tgt_homodimer.weight/2)
    ser_homodimer_r = get_score(df_combos["dg_r"], tgt_homodimer, wt=tgt_homodimer.weight/2)

elif str_kit == "taqman":
    ser_homodimer_f = get_score(df_combos["dg_f"], tgt_homodimer, wt=tgt_homodimer.weight/3)
    ser_homodimer_p = get_score(df_combos["dg_p"], tgt_homodimer, wt=tgt_homodimer.weight/3)
    ser_homodimer_r = get_score(df_combos["dg_r"], tgt_homodimer, wt=tgt_homodimer.weight/3)
ser_heterodimer = get_score(df_combos["dg_het"], tgt_heterodimer)

df_combos["score"] = ser_nprimers + ser_tmdiff + ser_homodimer_f + ser_homodimer_r + ser_heterodimer
df_combos.sort_values(by="score", inplace=True)

# ---------------------------------[ Get sequences from ambig fasta ]---------------------------------

if str_primerfile != "none":
    print("[rp] Retrieving sequences from ambig primer file")
    logger.info("Retrieving sequences from file: {0}".format(pth_primers))
    primers = {record.id: record.seq for record in SeqIO.parse(pth_primers, "fasta")}

    if str_kit == "endpoint": 
        # get primer sequences from dictionary, then reorder columns
        df_combos["seq_f"] = df_combos["id_f"].map(primers)
        df_combos["seq_r"] = df_combos["id_r"].map(primers)
        df_combos.insert(1, 'seq_f', df_combos.pop('seq_f'))
        df_combos.insert(5, 'seq_r', df_combos.pop('seq_r'))
    
    elif str_kit == "taqman":
        df_combos["seq_f"] = df_combos["id_f"].map(primers)
        df_combos["seq_r"] = df_combos["id_r"].map(primers)
        df_combos["seq_p"] = df_combos["id_p"].map(primers)
        df_combos.insert(1, 'seq_f', df_combos.pop('seq_f'))
        df_combos.insert(5, 'seq_p', df_combos.pop('seq_p'))
        df_combos.insert(9, 'seq_r', df_combos.pop('seq_r'))

# ===========================================[ RUN SUMMARY ]==========================================

#output file
print("[rp] Writing sorted psets to: {0}".format(pth_output.name))
logger.info("Writing to output file: {0}\n".format(pth_output))
df_combos.to_csv(pth_output, index=False)

logger.info("Program complete.")
print("[rp] Program complete.")