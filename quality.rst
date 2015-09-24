Short read quality and trimming
===============================

.. note::

   Reminder: if you're on Windows, you should install `mobaxterm <http://mobaxterm.mobatek.net/download.html>`__.

OK, you should now be logged into your Amazon computer! How exciting!

Prepping the computer
---------------------

Before we do anything else, we need to set up a place to work and
install a few things.

First, let's set up a place to work::

   sudo chmod a+rwxt /mnt

This makes '/mnt' a place where we can put data and working files.

----

Next, let's install a few things::

   sudo apt-get update
   sudo apt-get install -y trimmomatic fastqc python-pip python-dev

These are the Trimmomatic and FastQC programs, which we'll use below,
along with some software prerequisites that we'll need for other things
below.

Data source
-----------

We're going to be using a subset of data from the E. coli reference
published in `Chitsaz et al., 2011
<http://www.ncbi.nlm.nih.gov/pubmed/21926975>`__.  The data was
originally downloaded from
http://bix.ucsd.edu/projects/singlecell/nbt_data.html, E. coli
reference lane 6.

1. Copying in some data to work with.
-------------------------------------

We've loaded subsets of the data onto an Amazon location for you, to
make everything faster for today's work.  We're going to put the
files on your computer locally under the directory /mnt/data::

   mkdir /mnt/data

Next, let's grab the data set::

   cd /mnt/data
   curl -O -L https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/microbial-2015-09-24/ECOLI_R1.fastq.gz
   curl -O -L https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/microbial-2015-09-24/ECOLI_R2.fastq.gz

Now if you type::

   ls -l

you should see something like::

  -rw-rw-r-- 1 ubuntu ubuntu 418068323 Sep 24 15:04 ECOLI_R1.fastq.gz
  -rw-rw-r-- 1 ubuntu ubuntu 429978135 Sep 24 15:05 ECOLI_R2.fastq.gz

These are each 5m read subsets of the original data.  This is analogous
to what you would see if you did a HiSeq or MiSeq run on a bacterial
sample.

One problem with these files is that they are writeable - by default, UNIX
makes things writeable by the file owner.  Let's fix that before we go
on any further::

   chmod u-w *

We'll talk about what these files are below.

2. Copying data into a working location
---------------------------------------

First, make a working directory; this will be a place where you can futz
around with a copy of the data without messing up your primary data::

   mkdir /mnt/work
   cd /mnt/work

Now, make a "virtual copy" of the data in your working directory by
linking it in -- ::

   ln -fs /mnt/data/* .

These are FASTQ files -- let's take a look at them::

   less ECOLI_R1.fastq.gz

(use the spacebar to scroll down, and type 'q' to exit 'less')

Question:

* why are there R1 and R2 in the file names?

Links:

* `FASTQ Format <http://en.wikipedia.org/wiki/FASTQ_format>`__

3. FastQC
---------

We're going to use `FastQC
<http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__ to
summarize the data. We already installed 'fastqc' on our computer -
that's what the 'apt-get install' did, above.

Now, run FastQC on two files::

   fastqc ECOLI_R1.fastq.gz
   fastqc ECOLI_R1.fastq.gz

Now type 'ls'::

   ls -d *fastqc*

to list the files, and you should see::

   ECOLI_R1

We are *not* going to show you how to look at these files right now -
you need to copy them to your local computer to do that.  We'll show
you that tomorrow.  But! we can show you what they look like, because
I've made copiesd of them for you:

* `0Hour_ATCACG_L002_R1_001.extract_fastqc/fastqc_report.html <http://2015-may-nonmodel.readthedocs.org/en/latest/_static/0Hour_ATCACG_L002_R1_001.extract_fastqc/fastqc_report.html>`__
* `0Hour_ATCACG_L002_R2_001.extract_fastqc/fastqc_report.html <http://2015-may-nonmodel.readthedocs.org/en/latest/_static/0Hour_ATCACG_L002_R2_001.extract_fastqc/fastqc_report.html>`__

Questions:

* What should you pay attention to in the FastQC report?
* Which is "better", R1 or R2? And why?

Links:

* `FastQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__
* `FastQC tutorial video <http://www.youtube.com/watch?v=bz93ReOv87Y>`__

See `slide 39 and onwards <http://angus.readthedocs.org/en/2015/_static/2015-lecture2-sequencing.pptx.pdf>`__ for what BAD FastQC reports look like!

4. Trimmomatic
--------------

Now we're going to do some trimming!  We'll be using
`Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`__, which
(as with fastqc) we've already installed via apt-get.

The first thing we'll need are the adapters to trim off::

  curl -O -L http://dib-training.ucdavis.edu.s3.amazonaws.com/mRNAseq-semi-2015-03-04/TruSeq2-PE.fa

Now, to run Trimmomatic::

   TrimmomaticPE ECOLI_R1.fastq.gz ECOLI_R2.fastq.gz \
        ECOLI_R1.qc.fq.gz s1_se ECOLI_R2.qc.fq.gz s2_se \
        ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
        LEADING:2 TRAILING:2 \                            
        SLIDINGWINDOW:4:2 \
        MINLEN:25

You should see output that looks like this::

   ...
   Input Read Pairs: 5000000 Both Surviving: 4991513 (99.83%) Forward Only Surviving: 7422 (0.15%) Reverse Only Surviving: 782 (0.02%) Dropped: 283 (0.01%)
   TrimmomaticPE: Completed successfully

Capture the newly orphaned sequences like so::

   cat s1_se s2_se | gzip > ECOLI_orphans.qc.fq.gz

Questions:

* How do you figure out what the parameters mean?
* How do you figure out what parameters to use?
* What adapters do you use?
* What version of Trimmomatic are we using here? (And FastQC?)
* Do you think parameters are different for RNAseq and genomic data sets?
* What's with these annoyingly long and complicated filenames?
* why are we running R1 and R2 together?

Links:

* `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`__

5. FastQC again
---------------

Run FastQC again on the trimmed files::

   fastqc ECOLI_R1.qc.fq.gz
   fastqc ECOLI_R2.qc.fq.gz
   fastqc ECOLI_orphans.qc.fq.gz

And now view my copies of these files: 

* `0Hour_ATCACG_L002_R1_001.qc.fq_fastqc/fastqc_report.html <http://2015-may-nonmodel.readthedocs.org/en/latest/_static/0Hour_ATCACG_L002_R1_001.qc.fq_fastqc/fastqc_report.html>`__
* `0Hour_ATCACG_L002_R2_001.qc.fq_fastqc/fastqc_report.html <http://2015-may-nonmodel.readthedocs.org/en/latest/_static/0Hour_ATCACG_L002_R2_001.qc.fq_fastqc/fastqc_report.html>`__

Let's take a look at the output files::

   less ECOLI_R1.qc.fq.gz

(again, use spacebar to scroll, 'q' to exit less).

Questions:

* is the quality trimmed data "better" than before?
* Does it matter that you still have adapters!?

6. Interleave the sequences
---------------------------

Next, we need to take these R1 and R2 sequences and convert them into
interleaved form, for the next step.  To do this, we'll use scripts
from the `khmer package <http://khmer.readthedocs.org>`__, which we
need to install::

   sudo pip install -U setuptools
   sudo pip install khmer==2.0

Now, interleave the reads::

   interleave-reads.py ECOLI_R1.qc.fq.gz ECOLI_R2.qc.fq.gz --gzip \
      -o ecoli_ref-5m-trim.pe.fq.gz

and rename the orphans::

   cp ECOLI_orphans.qc.fq.gz ecoli_ref-5mtrim.se.fq.gz

Done!

Next: :doc:`assembling-ecoli`

