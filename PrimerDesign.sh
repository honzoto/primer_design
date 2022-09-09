#!/bin/bash

# =======================================[ DECLARING VARIABLES ]======================================
# set up default parameters
#DIR_SCRIPT="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
DIR_PROGRAM="/mnt/tank/bench/scripts/primer_design"
${DIR_PROGRAM}/header.sh

MAX_MISMATCH=2
DG_CUTOFF="-9"

# creating help function
display_help() {
    # taken from https://stackoverflow.com/users/4307337/vincent-stans
    echo "Usage: $0 [option...] " >&2
    echo
    echo "  [REQUIRED ARGUMENTS]"
    echo "  -p/--project        : <dir> project name (relative to 'projects' directory"
    echo "  -a/--alignment      : <file> multiple sequence alignment FASTA filename"
    echo "  -db/--database      : <dir> name of folder to be used as negative dataset"
    echo "  -mm/--mismatch      : <int> maximum allowed ambiguous nucleotides in a primer"
    echo "  -dg/--dg_cutoff     : <int> minimum delta-G requirement for heterodimers"
    echo "  -ic/--control       : <dir> name of internal controls database for filtering"
    echo
    echo "  [CASE 1]"
    echo "  Provide reference genome and coordinates to determine CDS region in MSA"
    echo "  -rf/--ref_fasta     : <file> reference FASTA filename"
    echo "  -rc/--ref_coords    : <file> coordinates of CDS regions in reference in CSV format"
    echo
    echo "  [CASE 2]"
    echo "  Provide set coordinates of CDS regions in MSA"
    echo "  -mc/--msa_coords    : <file> coordinates of CDS regions in MSA in CSV format"
    echo

    exit 1
}

# -------------------------------------[ Collect User Arguments ]-------------------------------------

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in        
        # case 1 - for determining new coordinates from MSA
        -rc|--ref_coords)
        REF_COORDS="$2"
        shift # past argument
        ;;
        -rf|--ref_fasta)
        REF_FASTA="$2"
        shift # past argument
        ;;

        # case 2 -forcing coordinates on MSA
        -mc|--msa_coords)
        MSA_COORDS="$2"
        shift # past argument
        shift # past value
        ;;

        # required args
        -p|--project)
        PROJECT="$2"
        shift # past argument
        shift # past value
        ;;
        -a|--alignment)
        ALIGNMENT="$2"
        shift # past argument
        ;;
        -db|--database)
        DATABASE="$2"
        shift # past argument
        ;;

        # optional args
        -mm|--mismatch)
        MAX_MISMATCH="$2"
        shift # past argument
        ;;
        -dg|--dg_cutoff)
        DG_CUTOFF="$2"
        shift # past argument
        ;;
        -ic|--controls)
        DIR_IC="$2"
        shift # past argument
        ;;

        -h|--help)
        display_help
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;
    esac
done

#
YELLOW=`tput setaf 3`
RESET=`tput sgr0`

DIR_PROJECT="/mnt/tank/bench/projects/"${PROJECT}
DIR_DATABASE="${DIR_PROGRAM}/databases/${DATABASE}"

PROG_PF="${DIR_PROGRAM}/bin/primer_figurator.py"
PROG_BF="${DIR_PROGRAM}/bin/blast_filterer.py"
PROG_DE="${DIR_PROGRAM}/bin/deltag_estimator.py"
PROG_SL="${DIR_PROGRAM}/bin/shortlister.py"

cd ${DIR_PROJECT}
FILE_LOG="log_primerdesign.txt"

# ---------------------------------------------[ Header ]---------------------------------------------

echo -e "\n${YELLOW}[PrimerDesign.sh] Starting Primer selection Pipeline...${RESET}"

# header block - reset log file
echo -e "Starting PrimerSelection Pipeline...\n" > ${FILE_LOG}
echo -e "INITIAL PARAMETERS:" >> ${FILE_LOG}
echo -e "Project directory\t\t\t[-p]: ${PROJECT}" >> ${FILE_LOG}
echo -e "Start time: $(date +'%Y/%m/%d') $(date +'%T')\n\n" >> ${FILE_LOG}

# --------------------------------------------[ Workflow ]--------------------------------------------

# Run primer_figurator.py
echo -e "\n------------------------------------[ Primer Selection ]------------------------------------\n" >> ${FILE_LOG}
PF_OUTFILE="primers_${MAX_MISMATCH}mm.fasta"
if [ ! -z "${REF_COORDS}" ]; then # if the variable is not empty
    echo "${YELLOW}Running primer_figurator.py with reference FASTA${RESET}"
    echo -e "[cmd] ${PROG_PF} -p ${PROJECT} -rc ${REF_COORDS} -rf ${REF_FASTA} \n\
    -a ${ALIGNMENT} -m ${MAX_MISMATCH} -o ${PF_OUTFILE}" >> ${FILE_LOG}
    ${PROG_PF} -p ${PROJECT} -o ${PF_OUTFILE} \
    -rc ${REF_COORDS} -rf ${REF_FASTA} -a ${ALIGNMENT} -m ${MAX_MISMATCH}
else
    echo "${YELLOW}Running primer_figurator.py with set coordinates${RESET}"
    echo -e "[cmd] ${PROG_PF} -p ${PROJECT} -mc ${MSA_COORDS} \n\
    -a ${ALIGNMENT} -m ${MAX_MISMATCH} -o ${PF_OUTFILE}" >> ${FILE_LOG}
    ${PROG_PF} -p ${PROJECT} -o ${PF_OUTFILE} \
    -mc ${MSA_COORDS} -a ${ALIGNMENT} -m ${MAX_MISMATCH}
fi

if [ ! "${MAX_MISMATCH}" = "0" ]; then
    PF_OUTFILE="primers_${MAX_MISMATCH}mm_expanded.fasta"
fi

# Run blast_filterer.py - with search option
echo -e "\n-------------------------------[ Negative Database Filtering ]------------------------------\n" >> ${FILE_LOG}
BF_OUTFILE="primers_${MAX_MISMATCH}mm_filtered.fasta"
echo -e "\n${YELLOW}[PrimerDesign.sh] Running blast_filterer.py${RESET}"
echo -e "[cmd] ${PROG_BF} -p ${PROJECT} -i ${PF_OUTFILE} \n\
    -o ${BF_OUTFILE} -db ${DATABASE} -ic ${DIR_IC}" >> ${FILE_LOG}
${PROG_BF} -p ${PROJECT} -i ${PF_OUTFILE} -o ${BF_OUTFILE} -db ${DATABASE} -ic ${DIR_IC}

# Run deltag_estimator.py
echo -e "\n----------------------------[ Homodimer & Heterodimer Filtering ]---------------------------\n" >> ${FILE_LOG}
echo -e "\n${YELLOW}[PrimerDesign.sh] Running deltag_estimator.py${RESET}"
echo -e "[cmd] ${PROG_DE} -p ${PROJECT} -i ${BF_OUTFILE} -ic ${DIR_IC}" >> ${FILE_LOG}
${PROG_DE} -p ${PROJECT} -i ${BF_OUTFILE} -ic ${DIR_IC}

# Run shortlister.py
echo -e "\n------------------------------[ Shortlisting Sets of Primers ]------------------------------\n" >> ${FILE_LOG}
echo -e "\n${YELLOW}[PrimerDesign.sh] Running shortlister.py${RESET}"
echo -e "${PROG_SL} -p ${PROJECT} -dg ${DG_CUTOFF}" >> ${FILE_LOG}
${PROG_SL} -p ${PROJECT} -dg ${DG_CUTOFF}

# ---------------------------------------------[ Footer ]---------------------------------------------

echo -e "[$(date +'%Y/%m/%d') $(date +'%T')] Program complete." >> ${FILE_LOG}

echo "${YELLOW}Program complete.${RESET}"