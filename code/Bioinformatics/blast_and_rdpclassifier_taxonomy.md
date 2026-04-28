#we will use blastn and a COI specific RDP classifier for taxonomy assignments of our ASVs

export BLASTDB=/media/mcrg/511fdc64-e3c3-4db7-9668-f8982f7782f9/blastdb/

#As of April 2026, using a conda environment for blast 2.17+
conda activate blast
#this includes other tools like blastx and perl scripts to update the reference database

## run with 1 top hit and 5 top hits

blastn -db nt_euk -query dna-sequences.fasta -max_target_seqs 1 \
-out 12Sblast_results.tsv \
-evalue 0.01 \
-outfmt "6 qseqid sseqid pident evalue length ssciname sblastname scomname" \
-num_threads 20

blastn -db nt_euk -query dna-sequences.fasta \
-max_target_seqs 5 \
-out 12Sblast_5results.tsv \
-evalue 0.01 \
-outfmt "6 qseqid sseqid pident evalue length ssciname sblastname scomname" \
-num_threads 20

#RDP classifier - mainly for COI. Change path to rRNAClassifier.properties as needed
rdp_classifier -Xmx32g classify \
-t /mnt/nvme1n1p1/eDNA_Data/RDPClassifier/mydata/rRNAClassifier.properties \
-o rdp.output dna-sequences.fasta
