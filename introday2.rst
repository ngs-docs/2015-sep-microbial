.. _intrornaseqday2:

RNA-seq and Differential Gene Expression in Bacteria
====================================================

Today we have a few objectives we would like to cover:

#. Set up a project for your data analysis for reproducible and robust analysis - :doc:`projectsetup`

#. Quality control and trimming of short RNA-seq reads - :doc:`qctrim`

#. Workflows for reference-based and reference-free transcriptome analysis - :doc:`aligncount` and :doc:`refvsnoref`

#. Differential gene expression and some potential pitfalls - :doc:`diffexp`


Learning goals
--------------

* Familiarize yourself with creating documentation and data organization strategies

* Perform and assess read trimming and quality

* Analyze reference and reference-free transcriptome analysis through available workflows

* Understand underlying assumptions and potential issues with differential gene expression


At this point, go ahead and start up new m3.xlarge EC2 instances, as you did
yesterday (see: :doc:`amazon/index`).

Additional Resource - Basic Linux/Unix commands
-----------------------------------------------

To refresh your memory on some basic Linux/Unix commands, we will cover the basic commands necessary to:

**1.** Move through folders

**2.** List the contents of a folder

**3.** Make new folders

**4.** Rename files/folders

**5.** Delete files/folders

.. csv-table::
   :header: " ", "Command", "What it does...", "Examples"
   :widths: 2, 8, 10, 40

   "**1.**", "cd", "Change directory/folder", "**>** cd ~ (this changes to your home directory); **>** cd .. (this goes back one folder)"
   "**2.**", "ls", "List the contents of a folder", "**>** ls"
   "**3.**", "mkdir", "Make a new directory/folder", "**>** mkdir NewFolder (this will make a new folder called 'NewFolder' in your current directory)"
   "**4.**", "mv", "Rename or move a file from one name to another", "**>** mv file1 file2 (this will rename/move file1 to file2)"  
   "**5.**", "rm", "Remove a file (add the -r flag to remove a folder)", "**>** rm file1 (remove file1); **>** rm -r folder1 (remove folder1)" 


**Command reference sheet**

.. image:: ./figures/linuxcoms.jpg
	:align: center
	:alt: Linux/Unix command list
	
(Ref. sheet from: `http://files.fosswire.com/2007/08/fwunixref.pdf <http://files.fosswire.com/2007/08/fwunixref.pdf>`__)

----

Next: :doc:`projectsetup`
