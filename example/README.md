# How to use AlFinder

##1) Downloading

Download ALL documents from https://github.com/tbigot/alFinder and save them in the same folder (e.g. NGS_analysis).
To perform this step, you can use git:

    git clone https://github.com/tbigot/alFinder.git

or the download  link on this repository website.


##2) Test the program with the example

Try the demo first in the `example` folder of the previously created folder (e.g. `NGS_analysis`)


1. In the example folder the user can find an example based on the method described in the following article:
> Large-scale genotyping by next generation sequencing: how to overcome the challenges to reliably genotype individuals?
Ferrandiz-Rovira M, Bigot T, Allainé D, Callait-Cardinal M-P, Radwan J, Cohas A.

 * the file `data.fna` contains the output data from the sequencing machine
 * the file `tags.csv` contains the tags
 * the files `alleles_locus1.fas`, `alleles_locus2.fas`, `alleles_locus3.fas` and `alleles_locus4.fas` contain the previously described alleles for four loci
 * the file `correspondenceData.csv` contains the correspondence file
 * the `output_demo` folder contains outputs from the DEMO settings  
 * the file `settings.ini` is modified according to the method described in the article and is ready to be use
 * the file `postprocessing.py` is modified according to the method described in the article and is ready to be use
2. Open the terminal (if python is not installed in your computer please download and install it on https://www.python.org/)
3. Use the command `cd` or `C:` depending on the operating system of your computer (linux/mac and Windows) to go in the directory where ALL downloaded documents have been saved (e.g. cd /home/NGS_analysis/DEMO)
4. Call `alFinder.py` file in the terminal (type `python alFinder.py` in the terminal)
5. Call `postprocessing.py` file in the terminal (type `python postprocessing.py` in the terminal)
6. Close the terminal
7. You can find the results in the previously created folder (e.g. `NGS_analysis/DEMO`)
    * the file `data_result.csv` is the obtained result after typing `python alFinder.py` in the terminal
    * the file `data_filteredResults.csv` is the obtained result after typing `python postprocessing.py` in the terminal

    
##3) Use it with your own data
You run succesfully the demo, now try it on your own data
1. Put your own data in the previously created folder (e.g. NGS_analysis):
 * a file with the output data from the sequencing machine
 * a file with the tags
 * a file/files with the previously described alleles (if no previously described alleles: an empty file should be present for each locus)
 * if necessary: the correspondence file
2. Modify `settings.ini` according to your needs, save and close it (see example in the `settings.ini` of the DEMO folder)
3. Modifiy `postprocessing.py` according to your needs, save and close it (see example in the `postprocessing.py` of the DEMO folder)
4. Open the terminal
5. Use the command `cd` or `C:` depending on the operating system of your computer (linux/mac and Windows) to go in the directory where ALL downloaded documents have been saved (e.g. `cd /home/NGS_analysis`)
6. Call `alFinder.py` file in the terminal (type `python alFinder.py` in the terminal)
7. Call `postprocessing.py` file in the terminal (type `python postprocessing.py` in the terminal)
8. You can find the results in the previously created folder (e.g. `NGS_analysis`)



