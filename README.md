# primer_design

Dockerfile currently under development.
get_primers.py, filter_db.py, generate_psets.py can be used individually, however, the directory referrals need to be changed into a portable style.

### get_primers.py
Creates high-throughput list of primers based on common regions of an <MSA> FASTA file. The MSA needs to be created through a third-party program (e.g. muscle). Output primers definition lines are formatted as such: ">gene|start-end|primer length|ambiguous bases~number". If there are more than one ambiguous bases, the ambiguous version of the primer fasta is written to ambig_<outputname>.fasta

[CASE 1] We use the entire MSA (we don't care about genes, any region of the MSA is valid target area)

[CASE 2] Coordinates pertaining to reference genome: We can also use a region (or multiple regions) from the MSA specified using a <ref coords> csv file, the coordinates should pertain to an <unaligned reference> that was used in the MSA. When you download a sequence on NCBI, they tell you the gene coordinates. You enter that into a coordinates file for the genes that you're interested in targeting. However, once you perform an MSA, these coordinates can change. Therefore, you need to tell get_primers.py what the original coordinates are in the original reference so that it can determine for you where the aligned genes are in the new MSA. In this case, the program will tell you your new coordinates as an output. 

[CASE 3] Coordinates pertaining to MSA: If you already know the MSA coordinates, you can use the -mc flag to specify the name of the CSV file. The new aligned genes will be written into a new folder called cds_seqs

USAGE CASE 1: get_primers,py -p <projectname> -o <primers_3mm.fasta> -a <MSA> [options]
USAGE CASE 2: get_primers.py -p <projectname> -o <primers_3mm.fasta> -a <MSA> -rc <ref coords> -rf <unaligned reference> [options]
USAGE CASE 3: get_primers.py -p <projectname> -o <primers_3mm.fasta> -a <MSA> -mc <msa coords> [options]

### filter_db.py
After generating the primers <primers_3mm.fasta> from get_primers.py, we need to filter out the primers that can potentially align to "off-targets". We do that with the help of BLAST searches. We keep all the BLAST databases (using makeblastdb) in the /database folder, and here, we specify which databases (<negative db>) we consider as off-targets.

Alternatively, we can also specify which databases our primers DO need to align to <positive db>. For example, COVID primers need to align to at least 90% of all reference genomes, we can specify that using the -pd flag.

USAGE: 
filter_db.py -p <projectname> -i <primers_3mm.fasta> -o <primers_3mm_filtered.fasta> -nd <negative db>

### generate_psets.py
Using all the primers available in the pool, generate_psets will determine which primers can potentially be used to create an amplicon (Taqman or Endpoint). First, primers are assembled by coordinates that would yield an appropriately-sized amplicon, then dG/Tm values are evaluated to see if the primers would form an appropriate PCR reaction.

These are the default params, but should be adjusted depending on the size of the primer pool you have available:
- Endpoint amplicons should have a length of 200-300 bp, and Taqman would have 60-150 bp
- The probe should have a Tm about 4-8 deg. C above the average of the forward and reverse. 
- The forward/reverse primers should have Tm within 2 deg. C of each other. 
- All delta-G values should be above -9.0

### rank_psets.py
Calculates a "score" for combinations from generate_psets.py, and ranks them
