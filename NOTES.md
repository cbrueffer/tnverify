Features
========

- Pipeline integration mode
- p-value whether two samples are paired
- input file which lists pairings?
- Bind everything into one binary
- possibility to add files to existing analysis
- Switch for storing mpileups in memory or on disk?


Correctness
===========

- test combinations of RNA-seq, DNA-seq, Exome-Seq
- test samples for depth to make sure they're sufficient for variant calling?
  => maybe better to just document this?
- threshold for calling two samples paired?
- What to do if reference or regions used to generate VCF files differ?


Related Tools
=============

VerifyBamID
http://genome.sph.umich.edu/wiki/VerifyBamID

IdCheck
http://eqtl.rc.fas.harvard.edu/idcheck/

Assessing Matched Normal and Tumor Pairs in Next-Generation Sequencing Studies 
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0017810

SNPrelate dendrogram functionality
http://bioinformatics.oxfordjournals.org/content/28/24/3326.short
