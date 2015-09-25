========================================
Assembling E. coli sequences with SPAdes
========================================

The goal of this tutorial is to show you the basics of assembly using
`the SPAdes assembler <http://bioinf.spbau.ru/spades>`__.

We'll be using data from `Efficient de novo assembly of single-cell
bacterial genomes from short-read data sets, Chitsaz et al., 2011
<http://www.ncbi.nlm.nih.gov/pubmed/21926975>`__.

Packages to install
===================

Download and insatll the SPAdes assembler::

   cd ~/
   wget http://spades.bioinf.spbau.ru/release3.6.0/SPAdes-3.6.0-Linux.tar.gz
   tar xvfz SPAdes-3.6.0-Linux.tar.gz
   export PATH=$PATH:$HOME/SPAdes-3.6.0-Linux/bin
   echo 'export PATH=$PATH:$HOME/SPAdes-3.6.0-Linux/bin' >> ~/.bashrc

as well as `Quast <http://quast.bioinf.spbau.ru/manual.html>`__,
software for evaluating the assembly against the known reference: ::

   cd
   curl -L http://sourceforge.net/projects/quast/files/quast-3.0.tar.gz/download > quast-3.0.tar.gz
   tar xvf quast-3.0.tar.gz

Getting the data
================

Now, let's create a working directory::

   cd /mnt
   mkdir assembly
   cd assembly

Copy in the E. coli data that you trimmed (:doc:`quality`)::

   ln -fs /mnt/work/ecoli_ref-5m-trim.*.fq.gz .

Running an assembly
===================

Now, let's run an assembly::

   spades.py --12 ecoli_ref-5m-trim.pe.fq.gz -s ecoli_ref-5m-trim.se.fq.gz -o spades.d

This will take about 15 minutes; it should end with the following output::

   * Corrected reads are in /mnt/assembly/spades.d/corrected/
   * Assembled contigs are in /mnt/assembly/spades.d/contigs.fasta (contigs.fastg)
   * Assembled scaffolds are in /mnt/assembly/spades.d/scaffolds.fasta (scaffolds.fastg)

Looking at the assembly
=======================

Run QUAST::

   ~/quast-3.0/quast.py spades.d/scaffolds.fasta -o report

and then look at the report::

   less report/report.txt

You should see::

   All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

   Assembly                   scaffolds
   # contigs (>= 0 bp)        152
   # contigs (>= 1000 bp)     80
   Total length (>= 0 bp)     4571384
   Total length (>= 1000 bp)  4551778
   # contigs                  89
   Largest contig             285527
   Total length               4558170
   GC (%)                     50.74
   N50                        133088
   N75                        67332
   L50                        12
   L75                        23
   # N's per 100 kbp          0.00

What does this all mean?

Comparing and evaluating assemblies - QUAST
===========================================

Download the true reference genome::

   cd /mnt/assembly
   curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoliMG1655.fa.gz
   gunzip ecoliMG1655.fa.gz

and run QUAST again::

   ~/quast-3.0/quast.py -R ecoliMG1655.fa spades.d/scaffolds.fasta -o report

Note that here we're looking at *all* the assemblies we've generated.

Now look at the results::

   less report/report.txt

and now we have a lot more information!  What all is there?

Adding in Nanopore data and doing a hybrid assembly
===================================================

One challenge with short read data like Illumina is that if there are
repeats in the genome, they can't be unambiguously resolved with short
reads.  Enter `long reads
<http://angus.readthedocs.org/en/2015/_static/Torsten_Seemann_LRS.pdf>`__,
produced by PacBio and Nanopore sequencing.  How much do long reads
improve the assembly?

Let's download some trial Nanopore data provided by Nick Loman::

   cd /mnt/assembly
   curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/microbial-2015-09-24/FC20.wf1.9.2D.pass.fasta.gz

Let's take a quick look at these sequences and try BLASTing them at
`NCBI <http://blast.ncbi.nlm.nih.gov/Blast.cgi>`__ -- note, you'll
want to use blastn, and choose "somewhat similar sequences" at the
bottom.  You can also restrict the BLAST search to E. coli MG1655.

Grab part of a sequence with ``gunzip -c FC20* | head`` and paste it
into BLAST.  What do you see?

Now, let's try adding them into the assembly by running SPAdes with
the Nanopore command line flag::

   spades.py --sc --12 ecoli_ref-5m-trim.pe.fq.gz -s ecoli_ref-5m-trim.se.fq.gz --nanopore FC20.wf1.9.2D.pass.fasta.gz -o nanopore-ecoli-sc

How'd we do? ::

   ~/quast-3.0/quast.py -R ecoliMG1655.fa nanopore-ecoli-sc/scaffolds.fasta -o n_report

Reference-free comparison
=========================

Above, we've been using the genome reference to do assembly
comparisons -- but often you don't have one. What do you do to
evaluate and compare assemblies without a reference?

One interesting trick is to just run QUAST with one assembly as a reference,
and the other N assemblies against it.  My only suggestion is to first
eliminate short, fragmented contigs from the assembly you're going to use
as a reference.

Let's try that, using ``extract-long-sequences.py`` from `khmer
<http://khmer.readthedocs.org>`__::

   extract-long-sequences.py -l 1000 nanopore-ecoli-sc/scaffolds.fasta > spades-long.fa

and then re-run QUAST and put the output in ``report-noref/report.txt``::

   ~/quast-3.0/quast.py -R spades-long.fa spades.d/scaffolds.fasta \
            nanopore-ecoli-sc/scaffolds.fasta -o report-noref

When you look at the report, ::

   less report-noref/report.txt

take particular note of the following -- ::

   Assembly                     spades.d_scaffolds  nanopore-ecoli-sc_scaffolds
   # contigs (>= 0 bp)          152                 15
   # contigs (>= 1000 bp)       80                  7
   Total length (>= 0 bp)       4571384             4643870
   Total length (>= 1000 bp)    4551778             4642289
   # contigs                    89                  7
   Largest contig               285527              3076878
   Total length                 4558170             4642289
   Reference length             4642289             4642289
      ...
   Misassembled contigs length  134677              0
   # local misassemblies        6                   0
   ...
   Genome fraction (%)          98.161              99.923
   Duplication ratio            1.000               1.001
   # mismatches per 100 kbp     3.36                0.00
