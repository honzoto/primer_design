#!/usr/bin/python3

import os, sys
from pathlib import Path
import numpy as np
from Bio import SeqIO

"""
PURPOSE:
This program creates a consensus sequence from an MSA

UPDATES:
2022-06-03: provided support for ambiguous bases

"""

str_program = Path(__file__).name
print("This is honzo's {0} to obtain consensus sequence from MSA.".format(str_program))
print("Working in directory: {0}".format(os.getcwd()))

# Default parameters
tf_degenerate = False
int_threshold = 1 # if less than int_threshold sequences have a mutation at a certain position, ignore it
str_directory = None

help_menu = """
----------------------------------[ HELP MENU ]---------------------------------

    USAGE: 
    python {p} -i <filename>

    PARAMETERS:
    -h/--help       : shows this menu

    User can select between "-i" for one input or "-d" for multiple inputs
    -i/--input      : <path> indicate single input MSA (.fasta) file
    -t/--threshold	: [default={t}] must have at least this number of bases
						from all records at a given position, otherwise, 
						the variation is ignored.
	-d/--degenerate	: [default={d}] use degenerate bases (IUPAC codes) for
						bases with multiple variants

--------------------------------------------------------------------------------
""".format(p=str_program, t=int_threshold, d=tf_degenerate)

for i, arg in enumerate(sys.argv):
	if arg == ("-i" or "--input"):
		pth_input = Path(sys.argv[i+1])
	elif arg == ("-t" or "--threshold"):
		int_threshold = int(sys.argv[i+1])
	elif arg == ("-d" or "--degenerate"):
		tf_degenerate = True

	elif arg == ("-h" or "--help"):
		print(help_menu); quit()

dic_ambs_rev = {"A":"A", "C":"C", "T":"T", "G":"G",
    "AC":"M", "AG":"R", "AT":"W", "CG":"S", "CT":"Y", "GT":"K",
    "ACG":"V", "ACT":"H", "AGT":"D", "CGT":"B", "ACGT":"N"}
dic_ambs = {"-":["-"], "A":["A"], "C":["C"], "T":["T"], "G":["G"],
	"M":["A","C"], "R":["A","G"], "W":["A","T"], "S":["C","G"],
    "Y":["C","T"], "K":["G","T"], "V":["A","C","G"], "H":["A","C","T"],
    "D":["A","G","T"], "B":["C","G","T"], "N":["A","C","G","T"]}

def get_mode(l: list):
	# return value with highest counts in a list
    return max(set(l), key = l.count)

def get_consensus(infile, threshold=1):
	print("Getting consensus sequence for {0}".format(infile))
	
	str_consensus = ""
	with open(infile, "r") as handle:
		# read each record, make it capitals, convert to a numpy array
		# when we have all the records, combine all the records into a 2D np array
		# we transpose the array so that we can read all bases at a given position at once
		reading_frame = np.array([list(str(r.seq).upper()) for r in SeqIO.parse(handle, "fasta")]).T

	for frame in reading_frame:
		expanded_frame = []
		for nuc in frame:
			expanded_frame += dic_ambs[nuc.upper()]

		if tf_degenerate:
			# we make two identical sets because we can't dynamically remove items from a set
			set_frame = set(expanded_frame)
			set_nucs = set(expanded_frame)
			if threshold > 0:
				# remove items from the set if it does not appear more than the threshold
				for nuc in set_frame:
					if list(frame).count(nuc) <= int_threshold:
						set_nucs.remove(nuc)
			lst_nucs = sorted(list(set_nucs))
		else:
			# get most common-occuring base
			lst_nucs = [get_mode(expanded_frame)]

		if lst_nucs == ["-"]:
			# if all nucleotides at that position is "-", then that frame doesn't have to be there
			continue
		elif "-" in lst_nucs or "N" in lst_nucs:
			str_consensus += "N"
		else:
			str_consensus += dic_ambs_rev["".join(lst_nucs)]
	return str_consensus

consensus = get_consensus(pth_input, int_threshold)
pth_output = pth_input.stem + "_cons" + ".fasta"
description = pth_input.stem + "_cons"

print("Writing to output file: {0}".format(pth_output))
with open(pth_output, "wt") as wf_out:
	wf_out.write(">{0}\n".format(description))
	wf_out.write(consensus+"\n")

print("Program complete.")
