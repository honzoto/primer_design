# -*- coding: utf-8 -*-


# Default python libraries
from pathlib import Path
import itertools

# Checking for dependencies
if 1 < 0:
    import pip
    lst_imports = ['PyQt5', 'pandas', 'biopython']
    for item in lst_imports:
        try: __import__(item); print(item+" already installed.")
        except ImportError: pip.main(['install', item])

from PyQt5 import QtCore, QtGui, QtWidgets
import tkinter as tk
from tkinter import filedialog

import pandas as pd
from Bio import SeqIO

str_version = "3.1"
"""
The purpose of this program is to estimate delta-G in a fasta file, generating a matrix of delta-g
values for every primer against every primer. The file can be made into a heatmap on Excel

version 1
- created in-house script that expands on all primers with ambiguous bases (dg_est.py)

version 2
- started development of GUI: main.py, created executable .pyc file
- moved main.py as additional program as part of PrimerDesign pipeline

"""
# GLOBAL VARIABLES
global dic_ambs, complement, dic_kmerdg, kmer_attr

class Kmer:
    def __init__(self, kmer, dg=-1):
        self.kmer = kmer
        self.dg = dg

    def set_tm(self, dh, ds, dg):
        self.tm_dh = dh
        self.tm_ds = ds
        self.tm_dg = dg

pth_kmerdg = Path("docs") / "experiment_2mers.csv"
pth_kmertm = Path("docs") / "santalucia-1996.csv"

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

complement = {"A":"T", "T":"A", "C":"G", "G":"C", "-":"N"}
dic_ambs = {"M":["A","C"], "R":["A","G"], "W":["A","T"], "S":["C","G"], \
        "Y":["C","T"], "K":["G","T"], "V":["A","C","G"], "H":["A","C","T"], \
        "D":["A","G","T"], "B":["C","G","T"], "N":["A","C","G","T"]}

root = tk.Tk()
root.withdraw()
root.call('wm', 'attributes', '.', '-topmost', True)

# ============================================[ FUNCTIONS ]===========================================

def revc(seq):
    s = seq.upper().replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c")
    return s[::-1].upper()

def count_kmers(seq, k=2):
    dic_counts = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in dic_counts:
            dic_counts[kmer] += 1
        else:
            dic_counts[kmer] = 1
    return dic_counts


def get_combinations(str_fullseq: str):
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
        print(i, str_subseq)
        lst_subseqs = [str_subseq[:-1]+nuc for nuc in dic_ambs[str_subseq[-1]]]
        lst_segments.append(lst_subseqs)
        j = i+1

    str_endseq = str_fullseq[j:]
    #print("e", str_endseq)

    lst_seqsegments = [seq for seq in itertools.product(*lst_segments)]
    lst_seqs = ["".join(seqsegment)+str_endseq for seqsegment in lst_seqsegments]
    return lst_seqs

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

## GUI
class Ui_frm_main(object):
    def setupUi(self, frm_main):
        frm_main.setObjectName("frm_main")
        frm_main.resize(540, 620)
        font = QtGui.QFont()
        font.setPointSize(10)
        frm_main.setFont(font)
        self.btn_execute = QtWidgets.QPushButton(frm_main)
        self.btn_execute.setGeometry(QtCore.QRect(280, 570, 120, 30))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.btn_execute.setFont(font)
        self.btn_execute.setObjectName("btn_execute")
        self.btn_quit = QtWidgets.QPushButton(frm_main)
        self.btn_quit.setGeometry(QtCore.QRect(410, 570, 120, 30))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.btn_quit.setFont(font)
        self.btn_quit.setObjectName("btn_quit")
        self.grb_files = QtWidgets.QGroupBox(frm_main)
        self.grb_files.setGeometry(QtCore.QRect(10, 10, 520, 111))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.grb_files.setFont(font)
        self.grb_files.setObjectName("grb_files")
        self.lbl_browseprimer = QtWidgets.QLabel(self.grb_files)
        self.lbl_browseprimer.setGeometry(QtCore.QRect(20, 30, 161, 30))
        self.lbl_browseprimer.setObjectName("lbl_browseprimer")
        self.btn_browseprimer = QtWidgets.QPushButton(self.grb_files)
        self.btn_browseprimer.setGeometry(QtCore.QRect(389, 32, 120, 26))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.btn_browseprimer.setFont(font)
        self.btn_browseprimer.setObjectName("btn_browseprimer")
        self.txt_primerfasta = QtWidgets.QLineEdit(self.grb_files)
        self.txt_primerfasta.setGeometry(QtCore.QRect(150, 30, 361, 30))
        self.txt_primerfasta.setObjectName("txt_primerfasta")
        self.lbl_profile = QtWidgets.QLabel(self.grb_files)
        self.lbl_profile.setGeometry(QtCore.QRect(20, 70, 161, 30))
        self.lbl_profile.setObjectName("lbl_profile")
        self.txt_profile = QtWidgets.QLineEdit(self.grb_files)
        self.txt_profile.setEnabled(False)
        self.txt_profile.setGeometry(QtCore.QRect(150, 70, 361, 30))
        self.txt_profile.setObjectName("txt_profile")
        self.txt_primerfasta.raise_()
        self.lbl_browseprimer.raise_()
        self.btn_browseprimer.raise_()
        self.lbl_profile.raise_()
        self.txt_profile.raise_()
        self.grb_params = QtWidgets.QGroupBox(frm_main)
        self.grb_params.setGeometry(QtCore.QRect(10, 130, 521, 281))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.grb_params.setFont(font)
        self.grb_params.setObjectName("grb_params")
        self.lbl_paramsets = QtWidgets.QLabel(self.grb_params)
        self.lbl_paramsets.setGeometry(QtCore.QRect(20, 30, 141, 25))
        self.lbl_paramsets.setObjectName("lbl_paramsets")
        self.lbl_targettype = QtWidgets.QLabel(self.grb_params)
        self.lbl_targettype.setGeometry(QtCore.QRect(20, 60, 141, 25))
        self.lbl_targettype.setObjectName("lbl_targettype")
        self.lbl_oligoconc = QtWidgets.QLabel(self.grb_params)
        self.lbl_oligoconc.setGeometry(QtCore.QRect(20, 100, 180, 25))
        self.lbl_oligoconc.setObjectName("lbl_oligoconc")
        self.lbl_naconc = QtWidgets.QLabel(self.grb_params)
        self.lbl_naconc.setGeometry(QtCore.QRect(20, 130, 180, 25))
        self.lbl_naconc.setObjectName("lbl_naconc")
        self.lbl_mgconc = QtWidgets.QLabel(self.grb_params)
        self.lbl_mgconc.setGeometry(QtCore.QRect(20, 160, 180, 25))
        self.lbl_mgconc.setObjectName("lbl_mgconc")
        self.lbl_dntpconc = QtWidgets.QLabel(self.grb_params)
        self.lbl_dntpconc.setGeometry(QtCore.QRect(20, 190, 180, 25))
        self.lbl_dntpconc.setObjectName("lbl_dntpconc")
        self.lbl_esttype = QtWidgets.QLabel(self.grb_params)
        self.lbl_esttype.setGeometry(QtCore.QRect(20, 230, 180, 25))
        self.lbl_esttype.setObjectName("lbl_esttype")
        self.txt_oligoconc = QtWidgets.QLineEdit(self.grb_params)
        self.txt_oligoconc.setGeometry(QtCore.QRect(240, 100, 113, 25))
        self.txt_oligoconc.setObjectName("txt_oligoconc")
        self.txt_naconc = QtWidgets.QLineEdit(self.grb_params)
        self.txt_naconc.setGeometry(QtCore.QRect(240, 130, 113, 25))
        self.txt_naconc.setObjectName("txt_naconc")
        self.txt_dntpconc = QtWidgets.QLineEdit(self.grb_params)
        self.txt_dntpconc.setGeometry(QtCore.QRect(240, 190, 113, 25))
        self.txt_dntpconc.setObjectName("txt_dntpconc")
        self.txt_mgconc = QtWidgets.QLineEdit(self.grb_params)
        self.txt_mgconc.setGeometry(QtCore.QRect(240, 160, 113, 25))
        self.txt_mgconc.setObjectName("txt_mgconc")
        self.lbl_oligounits = QtWidgets.QLabel(self.grb_params)
        self.lbl_oligounits.setGeometry(QtCore.QRect(360, 100, 50, 25))
        self.lbl_oligounits.setObjectName("lbl_oligounits")
        self.lbl_naunits = QtWidgets.QLabel(self.grb_params)
        self.lbl_naunits.setGeometry(QtCore.QRect(360, 130, 50, 25))
        self.lbl_naunits.setObjectName("lbl_naunits")
        self.lbl_mgunits = QtWidgets.QLabel(self.grb_params)
        self.lbl_mgunits.setGeometry(QtCore.QRect(360, 160, 50, 25))
        self.lbl_mgunits.setObjectName("lbl_mgunits")
        self.lbl_dntpunits = QtWidgets.QLabel(self.grb_params)
        self.lbl_dntpunits.setGeometry(QtCore.QRect(360, 190, 50, 25))
        self.lbl_dntpunits.setObjectName("lbl_dntpunits")
        self.cmb_paramsets = QtWidgets.QComboBox(self.grb_params)
        self.cmb_paramsets.setGeometry(QtCore.QRect(240, 30, 150, 25))
        self.cmb_paramsets.setObjectName("cmb_paramsets")
        self.cmb_targettype = QtWidgets.QComboBox(self.grb_params)
        self.cmb_targettype.setGeometry(QtCore.QRect(240, 60, 150, 25))
        self.cmb_targettype.setObjectName("cmb_targettype")
        self.cmb_esttype = QtWidgets.QComboBox(self.grb_params)
        self.cmb_esttype.setGeometry(QtCore.QRect(240, 230, 150, 25))
        self.cmb_esttype.setObjectName("cmb_esttype")
        self.txt_progress = QtWidgets.QPlainTextEdit(frm_main)
        self.txt_progress.setGeometry(QtCore.QRect(10, 420, 521, 121))
        self.txt_progress.setObjectName("txt_progress")

        self.retranslateUi(frm_main)
        self.set_connections()
        QtCore.QMetaObject.connectSlotsByName(frm_main)

        self.frm_main = frm_main

    def retranslateUi(self, frm_main):
        _translate = QtCore.QCoreApplication.translate
        frm_main.setWindowTitle(_translate("frm_main", "delta-G Estimator"))

        self.btn_execute.setText(_translate("frm_main", "E&xecute"))
        self.btn_quit.setText(_translate("frm_main", "&Quit"))
        self.grb_files.setTitle(_translate("frm_main", "Files for read"))
        self.lbl_browseprimer.setText(_translate("frm_main", "Primer fasta:"))
        self.btn_browseprimer.setText(_translate("frm_main", "Browse..."))
        self.lbl_profile.setText(_translate("frm_main", "Kmer profile:"))
        self.txt_profile.setText(_translate("frm_main", "/mnt/tank/bench/scripts/primer_design/docs/experiment_2mers.csv"))
        self.grb_params.setTitle(_translate("frm_main", "Input parameters"))
        self.lbl_paramsets.setText(_translate("frm_main", "Parameter sets:"))
        self.lbl_targettype.setText(_translate("frm_main", "Target type:"))
        self.lbl_oligoconc.setText(_translate("frm_main", "Oligo concentration:"))
        self.lbl_naconc.setText(_translate("frm_main", "Na+ concentration:"))
        self.lbl_mgconc.setText(_translate("frm_main", "Mg2+ concentration:"))
        self.lbl_dntpconc.setText(_translate("frm_main", "dNTP concentration:"))
        self.lbl_esttype.setText(_translate("frm_main", "Estimation algorithm:"))
        self.txt_oligoconc.setText(_translate("frm_main", "0.25"))
        self.txt_naconc.setText(_translate("frm_main", "20"))
        self.txt_dntpconc.setText(_translate("frm_main", "0.75"))
        self.txt_mgconc.setText(_translate("frm_main", "4.5"))
        self.lbl_oligounits.setText(_translate("frm_main", "uM"))
        self.lbl_naunits.setText(_translate("frm_main", "mM"))
        self.lbl_mgunits.setText(_translate("frm_main", "mM"))
        self.lbl_dntpunits.setText(_translate("frm_main", "mM"))
        self.txt_progress.setPlainText(_translate("frm_main", "--- welcome to Honzo\'s delta-G estmator application ---\n"))

        # Load items into combo-boxes
        self.cmb_targettype.addItems(['DNA', 'RNA'])
        self.cmb_targettype.setCurrentIndex(0)
        self.cmb_paramsets.addItems(['qPCR'])
        self.cmb_paramsets.setCurrentIndex(0)
        self.cmb_esttype.addItems(['single-base', '2-mer', '2-mer staggered'])
        self.cmb_esttype.setCurrentText('2-mer')

    def set_connections(self):
        self.btn_quit.clicked.connect(self.event_quit)
        self.btn_browseprimer.clicked.connect(self.event_selectfasta)
        self.btn_execute.clicked.connect(self.event_execute)

    def event_selectfasta(self):
        # Read all selected fasta files, separated by ';' character
        lst_primerfiles = filedialog.askopenfilenames(title="Select Primer FASTA files", filetypes=[("FASTA","*.fasta")])
        for i, str_primerfile in enumerate(lst_primerfiles):
            self.txt_progress.appendPlainText("File-{0}: {1}".format(i, str_primerfile))
        self.txt_primerfasta.setText(";".join(lst_primerfiles))

    def event_execute(self):
        # RUNNING MAIN FUNCTION
        self.btn_execute.setEnabled(False)

        # SETTING UP LOG FILE
        lst_primerfiles = [Path(f.strip().split(":")[-1]) for f in self.txt_primerfasta.text().split(";")]
        pth_logfile = lst_primerfiles[0].parent / "log_dg-main.txt"
        wf_log = open(pth_logfile, "wt")

        # Retrieving parameters
        wf_log.write("PARAMETERS\n")
        wf_log.write("K-mer delta-G file:   {0}\n".format(self.txt_profile.text()))
        wf_log.write("Parameter sets:       {0}\n".format(self.cmb_paramsets.currentText()))
        wf_log.write("Target type:          {0}\n".format(self.cmb_targettype.currentText()))
        wf_log.write("Oligo concentration:  {0} uM\n".format(self.txt_oligoconc.text()))
        wf_log.write("Na+ concentration     {0} mM\n".format(self.txt_naconc.text()))
        wf_log.write("Mg2+ concentration    {0} mM\n".format(self.txt_mgconc.text()))
        wf_log.write("dNTP concentration    {0} mM\n".format(self.txt_dntpconc.text()))
        wf_log.write("Est. algorithm        {0}\n\n".format(self.cmb_esttype.currentText()))

        # Getting K-mer profile from file
        pth_profile = Path(self.txt_profile.text())
        df_kmerdg = pd.read_csv(pth_profile)
        dic_kmerdg = dict(zip(df_kmerdg["kmer"], df_kmerdg["dg"]))

        # loop through the list first to write to log
        print("Getting list of primer files")
        for i, pth_fastafile in enumerate(lst_primerfiles):
            wf_log.write("File-{0}: {1}\n".format(i, pth_fastafile))
        
        for pth_fastafile in lst_primerfiles:
            wf_log.write("FASTA file: {0}\n".format(pth_fastafile.name))
            self.txt_progress.appendPlainText("\nReading FASTA file: {0}".format(pth_fastafile.name))

            # EXPANDING ON AMBIGUOUS PRIMERS
            dic_allseqs = {}
            with open(pth_fastafile) as handle:
                lst_records = list(SeqIO.parse(handle, "fasta"))
                for sio_raw in lst_records:
                    for i, str_expanded in enumerate(get_combinations(sio_raw.seq.upper()), start=1):
                        dic_allseqs[sio_raw.id+"-"+str(i)] = str_expanded

            wf_log.write("Created {0} non-ambiguous primers from {1} sequences.\n" \
                .format(len(lst_records), len(dic_allseqs)))

            # START ALIGNMENTS AND CALCULATIONS
            self.txt_progress.appendPlainText("Modelling and estimating Gibbs Free Energy")

            pth_outfile = pth_fastafile.parent / pth_fastafile.name.replace(".fasta", "_out.csv")
            pth_alignfile = pth_fastafile.parent / pth_fastafile.name.replace(".fasta", "_alignments.txt")
            with open(pth_alignfile, "wt") as wf_align, open(pth_outfile, "wt") as wf_out:
                wf_align.write("Alignment file: {0}\n".format(pth_alignfile))
                wf_align.write("Output matrix : {0}\n".format(pth_outfile))
                wf_out.write("primer,"+",".join(list(dic_allseqs.keys()))+"\n")

                for p1 in dic_allseqs:
                    wf_out.write("{0},".format(p1))
                    for p2 in dic_allseqs:

                        dic_dgresult = get_alignments(dic_allseqs[p1], dic_allseqs[p2][::-1], type=self.cmb_esttype.currentText(), dic_kmerdg=dic_kmerdg)
                        wf_out.write("{0},".format(dic_dgresult["dg"]))

                        str_towrite = "Primer 1: {n}".format(n=p1)
                        str_towrite += "seq: {0}\n".format(dic_allseqs[p1])
                        str_towrite += "Primer 2: {n}".format(n=p2)
                        str_towrite += "seq: {0}\n".format(dic_allseqs[p2])
                        str_towrite += "Estimated delta-G: {d}\n".format(d=round(dic_dgresult["dg"],2))
                        str_towrite += dic_dgresult["printout"]+"\n"
                        wf_align.write(str(str_towrite)+"\n")

                    wf_out.write("\n")

        wf_log.write("\nTask complete.")
        wf_log.close()

        self.txt_progress.appendPlainText("\nFinished.")
        self.btn_execute.setEnabled(True)

    def event_quit(self):
        print("Quitting.")
        self.frm_main.close()

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    frm_main = QtWidgets.QDialog()
    ui = Ui_frm_main()
    ui.setupUi(frm_main)
    frm_main.show()
    sys.exit(app.exec_())
