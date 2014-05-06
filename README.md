fgwas
=====

fgwas is a command line tool for integrating functional genomic information into a genome-wide association study (GWAS). The basic setup is as follows: you have performed a GWAS or a meta-analysis of many GWAS, and have identified tens of loci that influence the disease or trait (our approach works best if there are at least ~20 independent loci with p-values less than 5e-8). We set out to address the following questions:

1. Are these associations enriched in particular types of regions of the genome (coding exons, DNAse hypersensitive sites, etc.)?
2. Can we use these enrichments (if they exist) to identify novel loci influencing the trait?

The underlying model is described in: Pickrell JK (2014) [Joint analysis of functional genomic data and genome-wide association studies of 18 human traits](http://biorxiv.org/content/early/2014/01/22/000752). bioRxiv 10.1101/000752

The annotations used in this paper are available [from GitHub](https://github.com/joepickrell/1000-genomes).

##Installation##

###Dependencies###
fgwas depends on:

- the [GNU Scientific Library](http://www.gnu.org/software/gsl/)

- the [Boost Libraries](http://www.boost.org)

###Quick Start###
The most up-to-date release is: version 0.3.3. See ["Releases"](https://github.com/joepickrell/fgwas/releases) above.
After downloading fgwas-0.3.3.tar.gz at the link above, run:

>tar -xvf fgwas-0.3.3.tar.gz

>cd fgwas-0.3.3

>./configure

>make

This will create an executable file called fgwas in the src directory. The most common compilation error is that the configure script cannot find Boost or GSL. You may have to tell the script explicitly where to find them. For example, on OS X using macports, installations go to the non-standard path /opt/local/lib. To configure in this case, replace the above configure step with:

>./configure LDFLAGS=-L/opt/local/lib

Example data is available in the test_data/ directory. To ensure that fgwas is working, run:

> ./src/fgwas -i test_data/test_LDL.fgwas_in.gz -w ens_coding_exon

A user guide is available here: [fgwas v0.3.x User Guide](https://github.com/joepickrell/fgwas/blob/master/man/fgwas_manual.pdf)

Previous versions are available from the [Google Code repository](https://code.google.com/p/gwas/).

