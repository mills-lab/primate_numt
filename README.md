# primate_numt
all scripts used for primate numts

Insertion_flank_read_extract.pl is a perl script to extract flanking reads around the Numt insertion breakpoint and also its mate pair. 
This can then be used for doing a de novo assembly of the Numt insertion.

```
Usage: Insertion_flank_read_extract.pl (File with Numts and samples they are found in bedFormat.txt) (directory with BAM files)
```

Where you first enter the (name of the script) followed by the (input file in bed format) (the input file has chromosome breakpoint (for any insertion) start coordinate and end coordinate and all the samples the insertion was discovered in saperated by a comma(,)).
Followed by the (directory path) where all the BAM files are stored. Please look at the example for input file <bonobo_polymorphic_numts.txt>

>Example of input file - bonobo_polymorphic_numts.txt

Hotspot.R is an R script that is able map hotspot and do a permuation analysis to see if the hotspots are significant. 

```
Usage: Open the script and change names of three files <Reference_events.txt> <Polymorphic_events.txt> <Human_karyotype.txt>. 
You can then run different parts of the script based on your requirements (commented in the script). 
```

>Example of input files  - 1) Reference_events.txt 2) Polymorphic_events.txt 3) Human_karyotype.txt

## Citation
* Dayama, Gargi, Weichen Zhou, Javier Prado-Martinez, Tomas Marques-Bonet, and Ryan E. Mills. 2020. [Characterization of Nuclear Mitochondrial Insertions in the Whole Genomes of Primates](https://www.biorxiv.org/content/10.1101/2020.02.24.963504v2.abstract),
bioRxiv. `https://doi.org/10.1101/2020.02.24.963504`

* Dayama, Gargi, Sarah B Emery, Jeffrey M Kidd, and Ryan E. Mills. 2014. [The genomic landscape of polymorphic human nuclear mitochondrial insertions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4227756/pdf/gku1038.pdf),
Nucleic Acids Research, 2014, gku1038, `https://doi.org/10.1093/nar/gku1038`

* Weichen Zhou, Sarah B Emery, Diane A Flasch, Yifan Wang, Kenneth Y Kwan, Jeffrey M Kidd, John V Moran, Ryan E Mills,
[Identification and characterization of occult human-specific LINE-1 insertions using long-read sequencing technology](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1173/5680708), 
Nucleic Acids Research, 2019, gkz1173, `https://doi.org/10.1093/nar/gkz1173`
