#we will use blastn and a COI specific RDP classifier for taxonomy assignments of our ASVs

export BLASTDB=/mnt/nvme1n1p1/Cody/blastdb

## run with 1 top hit and 5 top hits

blastn -db nt_euk -query dna-sequences.fasta -max_target_seqs 1 \
-out 12Sblast_results.tsv \
-evalue 0.01 \
-outfmt "6 qseqid sseqid pident evalue length ssciname sblastname scomname" \
-num_threads 20

blastn -db nt_euk -query dna-sequences.fasta \
-max_target_seqs 5 \
-out 16Sblast_5results.tsv \
-evalue 0.01 \
-outfmt "6 qseqid sseqid pident evalue length ssciname sblastname scomname" \
-num_threads 20

rdp_classifier -Xmx32g classify \
-t /mnt/nvme1n1p1/eDNA_Data/RDPClassifier/mydata/rRNAClassifier.properties \
-o rdp.output dna-sequences.fasta
