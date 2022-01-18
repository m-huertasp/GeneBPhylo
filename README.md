# Gene-B Phylo — Working with OrthoFinder Output to analyse duplication and de novo events made automatic!


## What does Gene-B_Phylo do?
Gene-B Phylo is a platform that enables for fast and accurate analysis of OrthoFinder results. It infers the number of **duplicated** and **de novo** proteins as well as the **duplication events** and **de novo events** from a few output files from OrthoFinder. It allows you to perform a normalized analysis (events/branch length).

References:

[Emms, D.M. and Kelly, S. **(2015)** OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. **Genome Biology** 16:157](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2)

[Emms, D.M. and Kelly, S. **(2019)** OrthoFinder:Phylogenetic orthology inference for comparative genomics. **Genome Biology** 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y) 

### Gene-B Phylo in python
In order to run Gene-B Phylo you will need a specific python version - python 3.7
If you are running this programme with python 3.8 you may not be able to complete the analysis.

- Python 3.7 (https://www.python.org/)

#### Installing Dependencies
To perform an analysis Gene-B Phylo requires some dependencies to be installed in the python version you are using:

1. [**BioPython**](https://biopython.org/) 

2. [**ete3**](http://etetoolkit.org/) - In [Mol Biol Evol](Jaime Huerta-Cepas, François Serra and Peer Bork. "ETE 3: Reconstruction, analysis and visualization of phylogenomic data."  Mol Biol Evol (2016) doi:10.1093/molbev/msw046)

3. [**pandas**](https://pandas.pydata.org/)

Users can refer to the installation notes provided with these packages for more detailed instructions.

## Running Gene-B Phylo
1. You will use a file containing the following paths needed. If that file is not in the folder, it will be created as Gene-B_Phylo_settings.txt. 

*If you do not give any file, Gene-B Phylo will show the "help" option. The programme needs you to give a file that does no exist in order to create Gene-B_Phylo_settings.txt.

2. This settings.txt file is a text file with 7 lines, ordered in a specific manner:
- **Specie_Name** -> The Reference_Specie without any added words. If Reference_Specie = S_cerevisiae_proteins, Specie_Name = S_cerevisiae

- **Reference_Specie** -> The name of the specie you are interested in, as it is named in the OrthoFinder files.
*It is essential though, that all the other species have the same "added" word: s_cerevisiae_proteins, s_paradoxus_proteins...

- **OrthoFinder_path** -> Here you have to write the path to the OrthoFinder results folder in wich you have those files you are interested in analyse. To run OrthoEvolution it is necessary to use directly the Results folder created by OrthoFinder.

- **Proteins_fasta** -> All the proteins used by OrthoFinder of the Reference Specie. Here you have to write the path to the .fa file of the reference specie.

- **Output** -> You can name the output files in order to now the information within them. In a section below you can see how the output files are named using this.

- **Support** -> You will accept those duplication events above this number

*In the Duplication.tsv, Support = Proportion of expected species for which both copies of the duplicated gene are present.

** Gene-B Phylo reads line per line from ":" to the end. All you have to do is write the path to each file in its corresponding line and after ":" and without spaces.

The following files from OrthoFinder are used in this analysis:
1) Gene_Duplication_Events -> Duplications.tsv
2) Orthogroups -> Orthogroups.tsv
3) Orthogroups -> Orthogroups_UnassignedGenes.tsv
4) Orthologues -> Orthologues_Reference_Specie -> Reference_Specie_*.tsv
5) Species_Tree -> Species_Tree_rooted.txt
6) Species_Tree -> Species_Tree_node_labels.txt
7) Reference_Specie.fa (Reference Species proteins)

- Additional Information
For additional information about the OrthoFinder output files used in this analysis:
[OrthoFinder](https://github.com/davidemms/OrthoFinder)

3. If the Gene-B_Phylo_settings.txt file is already created and filled in the proper way the Gene-B Phylo analysis will be automatic.

4. The programme will print the date and time and its progress. It does not take more than five minutes in a vertebrates dataset of 17 species.


#To Run Gene-B Phylo on the Example Data type

`python3 Gene-B_Phylo1.1.py Gene-B_Phylo_settings.txt`


## What Gene-B Phylo provides
A Gene-B Phylo run produces a set of files describing an ordered by specie list of duplicated and de novo proteins. It also produces two Overall (de novo and duplication) files. In each one you can find a global summary of all the data extracted. 


### Results Files
1. **Date_Output_De_novo_Analysis.tsv** is a tab separated text file. First row contains the name of the gene, second one the Orthogroup in which the OrthoFinder classified that gene and two more rows containing the Clade were that gene is found and its corresponding node.

2. **Date_Output_De_novo_Overall.tsv** is a tab separated text file with six rows. The first two rows show each Node and Clade. The third one contains all de novo genes of that specific specie and the fourth all de novo events extracted. Then you can see the branch length (obtained from the OrthoFinder Species Tree) and a row named Normalized events, which accounts for (De novo events/(Branch length*1000)).
 
3. **Date_Output_Duplication_Analysis.tsv** tab separated file containing the same information as De novo Analysis but focused in Duplications.

4. **Date_Output_Duplication_Overall.tsv** tab separated file containing the same information as De novo Overall but focused in Duplications.
   
Enjoy!
