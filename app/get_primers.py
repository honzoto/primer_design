#!/usr/bin/python3

import sys, re, os
import glob, shutil, itertools
from pathlib import Path

from Bio import SeqIO
import logging

str_version = "1.3.10d"
str_program = Path(__file__).name
os.chdir(os.path.dirname(os.path.realpath(__file__)))
pth_base = Path(os.getcwd())

print("Welcome to honzo's {0} v{1} for extracting primer sequences.".format(str_program, str_version))

"""
< PURPOSE> 
The idea behind this program is to generate a list of primers based on the common areas of the genome.
Additionally, a coordinates file can be passed through the arguments if the user only wishes to look at certain
regions. A reference file will need to be supplied based on how which sequence the coordinates pertain to.

< HAHN'S NOTES >

v1  - rf"<some string> {var_name}\n" converts anything inside {} to actual variables,
        but interprets \n as a literal string, whereas "<some string {v}\n".format(v=var_name)
        will write what we usually write when it comes to log files
    - we will output our CDS genes into separate fasta files in cds_seqs in projects folder
    - NC001526.1 is a negative dataset - we removed it from the MSA
    - included two runtypes, one for CDS region isolation, one for primer generation
v2  - algorithm will now not write out subset primers (i.e. primers that are a subset of a larger primer),
        and also cut the ends of primers if not in [ACTG], given the length does not fall below length 
        criteria we now have primers with 0,1, and 2mm when we use maximum 2mm input argument
    - we now have a more detailed up-to-date help menu assigned directly in the script
    - using the script that Roy/Enaam developed (ambiguity_converter.py), we now have the option to 
        expand our ambiguities to only contain ACTG nucleotides at the end of the run
    - implemented "average" gc-based filtering to ease number of primers for downstream calculations
    - we can specify name of output file using -o/--output option
    - added colours to progress bar, implemented central log file
v3  - March 2022 pipeline re-design
    - implemented logging module for logging instead of manual writefile
    - using in-house developed ambiguous primer expander instead of external software
    - changed ambiguity identifier from "_" character to "~"
    - if a reference has a non-ACTG base region, no primers will be made there
    - program can now read multi-fastas in addition to single-line fastas
    - hz401 added support for ambiguous (non ACTG) bases in references

TODO
- bugs:
    - we aren't picking the correct starting locations - when running get_locations.py, we can see
        that at an alignment position, there are bases before that are the same

- consider changing primer-finding algorithm, for example, a primer of length 30 with 3 mm might not be
    as good as a primer of length 28 with 2mm, should we try to pick out the primer with the lowest
    mm ratio? how would we do that anyway? - done
- sometimes, we have 0 ambiguous bases in the primer, but still see |1 in the defline, this is because
    of the int_threshold parameter. If there is less than int_threshold ambiguous bases at that
    position, it will be noted as a mismatch, but does not show in the primer name
"""

# ===================================[ SETTING DEFAULT PARAMETERS ]===================================

tf_verbose = False # print warning messages that are in loops (can be excessive)
lst_runtype = ["cds", "pd"]
str_forcecoords = None # filepath for MSA coordinates
int_maxmm = 3
int_threshold = 0 # if 'n' or less genomes in MSA has a variant, ignore
lst_primerlens = [19, 22]
lst_gclimits = [30, 70]
str_refcoords = ""; str_reference = ""

str_helpmenu = \
"""
----------------------------------[ HELP MENU ]---------------------------------

 -h/--help          : shows this menu
 -p/--project       : <directory> name of the project; all subsequent '<files>'
                       should be referenced relative to the project directory
 -o/--output        : <string> name of output file (fasta)
 -t/--type          : [defualt={t}] select [cds/pd/cds,pd] for runtype
 -v/--verbose       : [default=False] print detailed warning and debugging
                       messages to the terminal

 [cds] options
 -rc/--ref_coords   : <file> coordinates of the CDS regions belonging to a 
                        reference fasta indicated by the flag '-r/--reference'
 -rf/--ref_fasta    : <file> reference FASTA where the coordinates pertain to
 -a/--alignment     : <file> multiple sequence alignment file, including with 
                        the reference used
 -mc/--msa_coords   : [default={fc}] <file> supply coordinates relative to MSA 
                        instead of reference. (disregards --ref_coords option)

 [pd] options
 -x/--max_mm        : [default={mm}] maximum number of mismatches allowed 
                        for a primer (mismatches include ambiguous bases)

--------------------------------------------------------------------------------
""".format(fc=str_forcecoords, mm=int_maxmm, t=lst_runtype)

# -------------------------------------[ Getting User Arguments ]-------------------------------------


print("Retrieving user arguments...")
lst_args = sys.argv
for i, arg in enumerate(lst_args):
    if arg == ("-p" or "--project"):
        str_project = lst_args[i+1]
    elif arg == ("-o" or "--output"):
        str_outfile = lst_args[i+1]
    elif arg == ("-t" or "--type"):
        lst_runtype = lst_args[i+1].split(",")

    # arguments for isolating CDS regions
    elif arg == ("-rc" or "--ref_coords"):
        str_refcoords = lst_args[i+1]
    elif arg == ("-rf" or "--ref_fasta"):
        str_reference = lst_args[i+1]
    elif arg == ("-a" or "--alignment"):
        str_msa = lst_args[i+1]
    elif arg == ("-mc" or "--msa_coords"):
        str_forcecoords = lst_args[i+1]

    # arguments for finding primers
    elif arg == ("-x" or "--max_mm"):
        int_maxmm = int(lst_args[i+1])
    elif arg == ("-h" or "--help"):
        os.system("../header.sh")
        print(str_helpmenu)
        quit()
    elif arg[0] == "-":
        print("[WARNING] unrecognized flag: {0}".format(lst_args[i]))
    
## SETTING UP GLOBAL VARIABLES
global dic_ambs, dic_ambs_rev
dic_ambs = {"A":["A"], "C":["C"], "G":["G"], "T":["T"], "-":["-"],
    "M":["A","C"], "R":["A","G"], "W":["A","T"], "S":["C","G"],
    "Y":["C","T"], "K":["G","T"], "V":["A","C","G"], "H":["A","C","T"],
    "D":["A","G","T"], "B":["C","G","T"], "N":["A","C","G","T"]}
dic_ambs_rev = {"A":"A", "C":"C", "T":"T", "G":"G",
    "AC":"M", "AG":"R", "AT":"W", "CG":"S", "CT":"Y", "GT":"K",
    "ACG":"V", "ACT":"H", "AGT":"D", "CGT":"B", "ACGT":"N"}

if (str_refcoords == "" or str_reference == "") and str_forcecoords == None:
    lst_runtype.remove("cds")

dir_project = Path("/projects/") / str_project
pth_log = dir_project / "PrimerDesign.log"
pth_alignment = dir_project / str_msa

pth_regex = dir_project / "reference_regex.txt"
pth_outcoords = dir_project / "coords_alignment.csv"


# ======================================[ CLASSES AND FUNCTIONS ]=====================================
class Strain:
    # class not being used efficiently, since we only will ever have one instance of a strain
    def __init__(self, name, seq=""):
        self.name = name
        self.gene = {}
        self.full_seq = seq

    def add_gene(self, gene_name, gene_coords):
        #print("[gp] subsetting gene {0} from {1} to {2}".format(gene_name, gene_coords[0], gene_coords[1]))
        self.gene_seq = self.full_seq[gene_coords[0]:gene_coords[1]]
        self.gene[gene_name] = self.Gene(self.gene_seq, gene_coords)

    def get_regex(self, outfile=None):
        """
        outfile: print regex sequences to a text file, similar to fasta format (for informational purposes)

        EXPLANATION:
        if the Strain is a reference genome, we want to convert the gene sequences into regular expressions
        this is so that we can search the aligned genome for where the gene is (if we can find it at all)
        for example, ATCGATCG does not match GCATGCA[ATC-GAT--CG], but A-*T-*C-*G-*... will
        """
        if outfile: open(outfile, "wt").close() # clear the existing regex log
        for gene_name in self.gene:
            logger.info("Generating regex for gene: {0}".format(gene_name))
            str_regex = ""
            for nuc in str(self.gene[gene_name].gene_seq): str_regex += nuc+"-*"
            self.gene[gene_name].regex = str_regex

            if outfile:
                with open(outfile, "at") as af_outfile:
                    af_outfile.write("{g}|start={s}|end={e}\n" \
                        .format(g=gene_name, s=self.gene[gene_name].coords[0], e=self.gene[gene_name].coords[1]))
                    af_outfile.write(str_regex+"\n")
                
    class Gene:
        def __init__(self, gene_seq, gene_coords):
            self.coords = gene_coords
            self.gene_seq = gene_seq

def get_gc(seq):
    # get the GC content of a sequence (average GC if it has ambiguities)
    seq = str(seq).upper()
    int_gc = seq.count('G') + seq.count('C')
    dic_ambs_gc = {'M':0.5,'S':1.0,'Y':0.5,'K':0.5,'V':2/3,'H':1/3,'D':1/3,'B':2/3,'N':0.5}

    for i in range(len(seq)):
        if seq[i] in dic_ambs_gc.keys():
            int_gc += dic_ambs_gc[seq[i]]

    return round(int_gc / len(seq) * 100, 1)

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


# =========================================[ PIPELINE START ]=========================================

# Setting up the logger 
print("[gp] Setting up logger: {0}".format(pth_log))
open(pth_log, "wt").close() # reset the log file
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
log_formatter = logging.Formatter("%(asctime)s-%(levelname)s: %(message)s", datefmt="[%Y-%m-%d %H:%M:%S]")
file_handler = logging.FileHandler(pth_log)
file_handler.setFormatter(log_formatter)
logger.addHandler(file_handler)

str_log = "This is honzo's {0} v{1} for determining primer sequences by common regions.\n"
str_log += "\t[cmd] {0}\n".format(" ".join(sys.argv))
str_log += "Project name            : {0}\n".format(str_project)
str_log += ""


# ---------------------------------[ Isolating CDS regions from MSA ]---------------------------------
if "cds" in lst_runtype:
    print("\nSTARING PART I: Isolating CDS regions from multiple sequence alignment")
    logger.info("\nIsolating for selected regions into cds_seqs directory...")

    # open already-aligned file from 3rd party software
    with open(pth_alignment) as handle:
        sio_alignments = list(SeqIO.parse(handle, "fasta"))

    # Keeping new alignment genes in cds_seqs directory relative to project folder
    try: 
        os.mkdir(dir_project / "cds_seqs")
    except FileExistsError:
        logger.warning("cds_seqs folder already exists. Overwriting existing files.")
        shutil.rmtree(dir_project / "cds_seqs")
        os.mkdir(dir_project / "cds_seqs")

    print("hz265 str_forcecoords:", str_forcecoords)
    if str_forcecoords:
        # -----------------------[ Retrieving pre-made MSA coordinates ]------------------------------

        logger.info(rf"Using preset coordinates from file: {str_forcecoords}")
        dic_coordinates = {}
        dir_forcecoords = dir_project / str_forcecoords
        with open(dir_forcecoords, "rt") as wf_forcecoords:
            wf_forcecoords.readline() # skip header
            for line in wf_forcecoords:
                lst_line = line.strip().split(",")
                dic_coordinates[lst_line[0]] = [int(c) for c in lst_line[1:]]

    else:
        # -----------------[ Calculating MSA coordinates from reference file ]------------------------

        pth_coords = dir_project / str_refcoords
        pth_reference = dir_project / str_reference
        print(rf"Determining CDS regions using reference file: {pth_reference}")
        logger.info("Reading reference file: {0}".format(pth_reference))

        # assuming only one FASTA record in the reference file
        with open(pth_reference) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                print("[gp] Reading reference:", record.id)
                stn_reference = Strain(record.id, record.seq)

        logger.info("Reading coordinates file: {0}".format(pth_coords))
        with open(pth_coords, "rt") as rf_coords:
            rf_coords.readline() # skip header
            for geneline in rf_coords:
                lst_geneline = geneline.split(",")
                # subtract one because gene position 1 is actually index 0 - this is the only time we'll do this
                stn_reference.add_gene(lst_geneline[0], [int(lst_geneline[1])-1, int(lst_geneline[2])])

        logger.info("\nGenerating regex sequences for reference genes.")
        stn_reference.get_regex()

        # now, we need to find the reference in our alignment file and locate the new coordinates
        dic_coordinates = {} # dictionary of subsets
        print(rf"[gp] getting gapped coordinates from alignment file for name: {stn_reference.name}")
        for record in sio_alignments:
            if str(record.id).find(stn_reference.name) >= 0:
                print(rf"[gp] found reference ({stn_reference.name}) in aligment file")
                str_coordinfo = "Using reference {0} as alignment\n".format(record.id)
                for str_gene in stn_reference.gene:
                    str_regex = stn_reference.gene[str_gene].regex
                    re_match = re.search(rf"{str_regex}", str(record.seq))
                    dic_coordinates[str_gene] = [re_match.start(), re_match.end()]
                    str_coordinfo += "gene {0}: {1}\n".format(str_gene, re_match)
                logger.info(str_coordinfo)
                break # this loop is only to get the information for the reference

        # write alignment coordinates to file (for us to use it next time)
        print(rf"[gp] writing new MSA coordinates to file: {pth_outcoords}")
        with open(pth_outcoords, "wt") as wf_outcoords:
            wf_outcoords.write("gene,start,end\n")
            for str_gene in dic_coordinates:
                wf_outcoords.write("{0},{1},{2}\n" \
                    .format(str_gene, dic_coordinates[str_gene][0]+1, dic_coordinates[str_gene][1]))

    # ------------------------[ Using new coordinates file to subset MSA ]----------------------------
    logger.info("\nUsing regex-based coordinates to extract genes from MSA.")
    lst_cdsfiles = [] # list of all the subset gene files

    for gene in dic_coordinates:
        pth_genefile = dir_project / "cds_seqs" / (gene+".fasta")
        print(rf"[gp] Creating aligned fasta file for gene: {gene}")
        start = dic_coordinates[gene][0]
        end = dic_coordinates[gene][1]

        with open(pth_genefile, "wt") as wf_genefile:
            for record in sio_alignments:
                wf_genefile.write(">"+record.id+"\n")
                wf_genefile.write(str(record.seq)[start:end]+"\n")

        logger.info("Created MSA for gene: {0}".format(pth_genefile))
        lst_cdsfiles.append(pth_genefile)


# ======================[ SEARCHING FOR PRIMERS IN MULTIPLE-SEQUENCE ALIGNMENT ]======================

if "pd" in lst_runtype:
    if str_outfile == "": str_outfile = "primers_{0}mm.fasta".format(int_maxmm)
    str_outfile = "ambig_"+str_outfile

    print("\nSTARTING PART II: Generating primers and writing to file {0}".format(str_outfile))

    pth_outfile = dir_project / str_outfile
    lst_alignmentfiles = glob.glob(str(dir_project / "cds_seqs")+"/*.fasta")
    if len(lst_alignmentfiles) == 0:
        # meaning we did not indicate CDS regions
        try: os.mkdir(dir_project / "cds_seqs")
        except: print("cds_seqs folder already exists... files may be overwritten")
        str_msacopy = str(dir_project / "cds_seqs") + "/" + str_msa.split("_")[0]+".fasta"
        shutil.copyfile(dir_project / str_msa, str_msacopy)
        lst_alignmentfiles = [str_msacopy]

    logger.info("\nSearching for primers and writing to output: {0}".format(pth_outfile))
    int_primersmade = 0
    int_poorgc = 0

    with open(pth_outfile, "wt") as wf_outfile:
        for str_genealignment in lst_alignmentfiles:
            logger.info("\nGetting primers from file: {0}".format(str_genealignment.split("/")[-1]))
            str_gene = str_genealignment.split("/")[-1].split(".")[0]

            lst_seqs = []
            lst_ids = [] 

            print(rf"[gp] Reading alignment file: {Path(str_genealignment).name}")
            for record in SeqIO.parse(str_genealignment, "fasta"):
                lst_ids.append(record.id)
                lst_seqs.append(record.seq)

            # --------------------------[ Searching for primers ]-------------------------------------
            n = -1 # nucleotide frame counter
            int_maxrun = 0
            int_totalframes = len(lst_seqs[0]) - lst_primerlens[1] - 1 # subtract one because of \n at end
            lst_maxrunloc = [0]
            str_previousprimer = "" # to filter out subset primers
            str_previousdefline = "A0|1-2|1|0"

            while n < int_totalframes:
                n += 1
                int_mm = 0
                str_primer = ""

                pct_bestmm = 1.00 # start at 100% mismatch, find lowest in sequence
                str_bestprimer = ""

                for p in range(lst_primerlens[1]): # max desired length of primers
                    lst_nucs = [] # lst_nucs is a list of all base occurrences
                    for s in range(len(lst_seqs)):
                        # hz401 - added support for ambiguous bases in reference
                        base = lst_seqs[s][n+p]
                        lst_nucs += dic_ambs[base]

                    # convert lst_nucs to unique bases only
                    set_nucs = set([nuc for nuc in lst_nucs if lst_nucs.count(nuc) > int_threshold])
                    lst_nucs = sorted(list(set_nucs)) # sort cos we'll eventually join them and reference dictionary

                    if "-" in lst_nucs:
                        if len(str_primer) == 0: break # we don't want primers starting with an 'N'
                        str_primer += "N"
                        int_mm += 1
                    else:
                        try:
                            str_key = dic_ambs_rev["".join(lst_nucs)]
                            str_primer += str_key
                            if str_key not in ["A","C","T","G"]: 
                                int_mm += 1
                        except KeyError:
                            if tf_verbose:
                                logger.warning("Unrecognized base(s) at position {0}: {1}".format(n+p, lst_nucs))
                            str_primer += "N"
                            int_mm += 1
                    
                    if int_maxrun < p: 
                        int_maxrun = p
                        lst_maxrunloc = [n]
                    elif int_maxrun == p:
                        lst_maxrunloc.append(n)

                    pct_mm = int_mm / len(str_primer)
                    # using <= operator because the primer length will get longer while still having 0 mm
                    if pct_mm <= pct_bestmm and len(str_primer) >= lst_primerlens[0]:
                        str_bestprimer = str_primer
                        pct_bestmm = pct_mm

                    # we have two possible cases for the primer extension to terminate:
                    # that's if the number of mismatches is greater than the maximum allowed,
                    if int_mm > int_maxmm: break # here to treat 0-mm primers
                    elif int_mm >= int_maxmm and p >= lst_primerlens[0]:
                        if p < lst_primerlens[0]: break # primer length not long enough
                        str_defline = ">{g}|{s}-{e}|{l}|{x}" \
                            .format(g=str_gene, s=n+1, e=n+len(str_bestprimer), l=len(str_bestprimer), x=int_mm)

                        # we don't want our current primer to be a subset of the previous primer
                        if str_previousprimer.find(str_bestprimer.strip()) < 0:
                            flt_gc = get_gc(str_bestprimer)
                            if flt_gc < lst_gclimits[0] or flt_gc > lst_gclimits[1]:
                                int_poorgc += 1
                            else:
                                int_primersmade += 1
                                wf_outfile.write(str_defline+"\n")
                                wf_outfile.write(str_bestprimer+"\n")
                                str_previousprimer = str_primer

                            break # TODO - unsure to put this statement here, or 1 tab ahead

                    elif p == lst_primerlens[1] - 1:
                        # two consecutive primers of the same length should not pass the .find() function
                        # so we don't need to check if the new primer is a subset of the old
                        
                        try:
                            while str_bestprimer[-1] not in ["A","C","T","G"] and len(str_bestprimer) > lst_primerlens[0]: 
                                str_bestprimer = str_bestprimer[:-1]
                                int_mm -= 1

                            # reached full length primer"
                            flt_gc = get_gc(str_bestprimer)
                            if flt_gc < lst_gclimits[0] or flt_gc > lst_gclimits[1]:
                                int_poorgc += 1
                            else:
                                int_primersmade += 1
                                str_defline = ">{g}|{s}-{e}|{l}|{x}" \
                                    .format(g=str_gene, s=n+1, e=n+p+1, l=len(str_bestprimer), x=int_mm)
                                wf_outfile.write(str_defline+"\n")
                                wf_outfile.write(str_bestprimer+"\n")
                                str_previousprimer = str_bestprimer
                        except:
                            if not str_primer.strip() == "":
                                print("Cannot process primer: {0}".format(str_primer))
                        break

        
    logger.info("Gene {g}: Longest primer: {p}-nt".format(g=str_gene, p=int_maxrun+1))
    logger.info("{0} possible primers were eliminated due to poor GC content".format(int_poorgc))
    logger.info("{0} primers were written to file: {1}.".format(int_primersmade, pth_outfile.name))
    print("---> {0} primers were found with {1} mismatches or less.".format(int_primersmade, int_maxmm))

    # ======================[ COVERTING AMBIGUOUS SEQUENCES TO MULTIPLE FASTAS ]======================

    pth_expanded = pth_outfile.parent / str_outfile.replace("ambig_", "")

    if int_maxmm > 0:
        print("[gp] Converting ambiguous codes to ACTG sequences.")
        logger.info("\nConverting ambiguous codes to ACTG sequences.")

        with open(pth_outfile) as handle:
            sio_ambigs = list(SeqIO.parse(handle, "fasta"))

        int_expandedseqs = 0
        with open(pth_expanded, "wt") as wf_expanded:
            for record in sio_ambigs:
                lst_expanded = get_combinations(record.seq)
                for i, seq in enumerate(lst_expanded):
                    int_expandedseqs += 1
                    wf_expanded.write(">{0}~{1}\n".format(record.id, i+1))
                    wf_expanded.write("{0}\n".format(seq))

        logger.info("{0} expanded primers written to file: {1}".format(int_expandedseqs, pth_expanded.name))
        print("---> {0} non-ambiguous primers generated.".format(int_expandedseqs))

    else:
        os.rename(pth_outfile, pth_expanded)
        logger.info("\nRenaming {0} to {1}".format(pth_outfile.name, pth_expanded.name))

logger.info("\n"+"-"*80+"\n")
print("[gp] Program complete.")
