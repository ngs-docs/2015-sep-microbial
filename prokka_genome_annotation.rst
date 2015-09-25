================================================
Bacterial genome annotation using Prokka
================================================

After you have de novo assembled your genome sequencing reads into contigs,
it is useful to know what genomic features are on those contigs. The process
of identifying and labelling those features is called genome annotation.

In this tutorial you will:

1. Download and install Prokka
2. Annotate a FASTA file of contigs
3. Visualize the annotation using Artemis

The instructions below will work on a Linux server (eg. EC2 instance),
a Linux desktop, and even directly on your Mac OS X Laptop.

Download and install Prokka
===========================

Prokka is simple to install because it comes bundled with all its dependencies:

::

  cd
  sudo apt-get -y install git bioperl libxml-simple-perl ncbi-blast+
  git clone https://github.com/tseemann/prokka.git
  export PATH=$PWD/prokka/bin:$PATH
  prokka --setupdb
  prokka --version

Apply Prokka to your just-assembled E. coli genome
==============================================

Prokka is a pipeline script which coordinates a series of genome feature predictor tools and sequence similarity
tools to annotate the genome sequence (contigs).

::

  cd /mnt/assembly
  prokka --outdir anno --prefix prokka --centre X --compliant spades-long.fa

::

  cat ./anno/prokka.txt

How many genes did Prokka find in the contigs?


Install Artemis
===============

Artemis is a graphical Java program to browse annotated genomes.
It is a a bit like IGV but specifically designed for bacteria.
You will need to install this on your desktop computer.
You could run it remotely over SSH using X11 forwarding from Amazon
but it is probably too slow to be useful.

Download: https://www.sanger.ac.uk/resources/software/artemis/#download


Viewing the annotated genome
============================

* Start Artemis
* Go to "File" -> "SSH File Manager"
* Type in the IP number of your Amazon EC2 instance
* Browse to the "anno/prokka.gff" file
