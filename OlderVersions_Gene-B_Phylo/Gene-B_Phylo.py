###############################################################
## Marta Huertas
## Evolutionary Genomics Group - Mar Albà group
## Research Unit on Biomedical Informatics (GRIB)
## UPF/IMIM, Barcelona Biomedical Research Park (PRBB)
###############################################################

#IMPORTING PACKAGES
import csv
import glob
import pandas as pd
from Bio import Phylo
from Bio import SeqIO
import ete3
from ete3 import Tree
import os.path as path
import sys
from datetime import datetime
from datetime import date

#FUNCTIONS
#With this function I am trying to get rid of all those duplicated genes in a Clade. It also renames the gene names in orden to have them all in the same format.
def extract_duplication_events (Clade_duplication_table):
    Total_dupl_genes = []
    Genes_1 = []
    Genes_2 = []
    for gene in Clade_duplication_table['Genes 1']:
        Genes_1.append(gene)
    for gene in Clade_duplication_table['Genes 2']:
        Genes_2.append(gene)
    Genes_1 = ','.join(Genes_1)
    Genes_1 = Genes_1.replace(' ','')
    Genes_2 = ','.join(Genes_2)
    Genes_2 = Genes_2.replace(' ','')
    Unique_Genes = Genes_1.split(',') + Genes_2.split(',')
    for gene in Unique_Genes:
        if Reference_Specie in gene:
            gene = gene.replace(Reference_Specie, '')
            gene = gene[1:]
            if gene not in Total_dupl_genes:
                Total_dupl_genes.append(gene)
    return Total_dupl_genes

#This function is used to obtain the species names from the oldest one to our Reference_Specie
def get_subtree (Reference_Specie, Tree):
    Clades = []
    for i in range(0,2):
        if Reference_Specie not in Tree[i].get_leaf_names():
            Clades.append(Tree[i].get_leaf_names())
        else:
            Subtree = Tree[i].get_children()
    return Clades, Subtree

#The following function is used to shorten the names of the species so it is easier to read them.
def short_names (Species):
    Clades_short = []
    for species in Species:
        if '_proteins' in species:
            Clades_short.append(species.replace('_proteins', ''))
        if 'v2' in species:
            Clades_short.remove(species.replace('_proteins', ''))
            Clades_short.append(species.replace('_v2_proteins', ''))
        if 'v3' in species:
            Clades_short.remove(species.replace('_proteins', ''))
            Clades_short.append(species.replace('_v3_proteins', ''))
    return Clades_short

#The following function is used to transform a DataFrame into a dictionary(key=Gene: value=Orthogroup)
def orthogroups_todict (DataFrame, Species):
    Orthogroup_Specie = DataFrame[['Orthogroup',Species]].dropna()
    Orthogroups_DICT = Orthogroup_Specie.set_index('Orthogroup')[Species].to_dict()
    Orthogroup_Genes = {}
    for element in Orthogroups_DICT:
        Orthogroup_Genes[element] = Orthogroups_DICT[element].replace(' ','').split(',')
    Gene_Orthogroup = {}
    for element in Orthogroup_Genes:
        for gene in Orthogroup_Genes.get(element):
            Gene_Orthogroup[gene] = element
    return Gene_Orthogroup


print(sys.argv[0])
if (len(sys.argv) != 2) or (sys.argv[1] == '-h'):
    print('USAGE:')
    print('Run OrthoEvolution.py from OrthoFinder Outputs directories ordered in a settings file')
    print('This programme requires, then, a settings.txt file to work')
    print('If this file is not found it will be created')
    print('For the settings file to be created you will need to use a random file which does not exist in your OrthoEvolution run')
    print('More information is to be found in the README.md file')
    print('This programme works with python version 3.7')
    print('----------------------------------------------------------------------------------------------------')
    print('Four output files will be created:')
    print('Date-Output_Duplications_Overall.tsv, Date-Output_Duplication_Analysis.tsv, Date-Output_De_novo_Overall.tsv, Date-Output_De_novo_gene_Analysis.tsv')
    
    
else:    
    ### First we need to get the information from a "Settings_File". There we have all the paths of all the files we are going to need.
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), ' : Starting ', sys.argv[0])
    print('Cheking required settings.txt file exists') 
    File = sys.argv[1]
    File_paths = []
    if path.isfile(File) == True:
        print('File found - ok')
        with open (File, 'r') as file:
            for line in file.readlines():
                File_paths.append(line.rstrip().split(':')[1])
            file.close()

        #Name of the reference specie we are going to use
        Reference_Specie = File_paths[0]
        Reference_Specie_short = Reference_Specie.replace('_proteins','')
        #File called Duplications.tsv
        Duplications_File = File_paths[1]
        #Species labeled and non labeled trees
        Species_Tree = File_paths[2]
        Non_labeled_tree = File_paths[3]
        #Folder containing all the files of the Orthologues of the reference specie
        Orthologues_Folder = File_paths[4]
        #All the proteins of the reference specie used
        Proteins_Fasta = File_paths[5]
        #Files containing the orthogroups (Orthogroups.tsv) and the orthogroups of the Unassigned genes(Orthogroups_UnassignedGenes.tsv)
        Orthogroups_File = File_paths[6]
        Orthogroup_UnassignedGenes = File_paths[7]
        Output = File_paths[8]
        Support = float(File_paths[9])

        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), ' : Running Duplication Analysis')  
        
#-----------------------------------------------------------------------------------------------#
#DUPLICATIONS
#-----------------------------------------------------------------------------------------------#
        #DataFrame showing all the duplications provided by OrthoFinder.
        DUPLICATIONS_DF = pd.read_csv(Duplications_File, delimiter = "\t", header=0)
        
        #We are going to get the clades we are interested in (related to the reference specie) directly from the Species_Tree
        S_tree = Phylo.read(Species_Tree, "newick")

        Tree_Path = []
        for clade in S_tree.trace('N0', Reference_Specie):
            c = str(clade)
            Tree_Path.append(c)
        Tree_Path = list(reversed(Tree_Path))
        Tree_Path.remove(Reference_Specie)

        #Names from the species tree -> With this we obtain a list of lists (if some species are agrupated they will be inside a list)
        Tree = ete3.Tree(Non_labeled_tree)
        Clades = []
        for i in range(0,len(Tree_Path)):
            if i == 0:
                Clades += get_subtree(Reference_Specie, Tree.get_children())[0]
                Subtree = get_subtree(Reference_Specie, Tree.get_children())[1]
            else:
                Clades += get_subtree(Reference_Specie, Subtree)[0]
                Subtree = get_subtree(Reference_Specie, Subtree)[1]
        Clades=list(reversed(Clades))

        #It is useful though to have all the species names in just one list:        
        Clades_Only = []
        for element in Clades:
            if len(element) == 1:
                Clades_Only += element
            else:
                for clade in element:
                    Clades_Only.append(clade)
        Clades_Only = short_names(Clades_Only)

        Species_Path = []
        for element in Clades:
            if len(element) == 1:
                Species_Path += short_names(element)
            if len(element) > 1:
                Species_Path.append(element[1].replace('_proteins','') + '_group')

        #This code is used to now the node in wich the novo gene is originated
        Clade_Path = {Reference_Specie_short:Reference_Specie_short}
        for i in range(0,len(Tree_Path)):
            Clade_Path[Species_Path[i]] = Tree_Path[i]

        Tree_Path = list(reversed(Tree_Path))
        Tree_Path.append(Reference_Specie)
        Tree_Path = list(reversed(Tree_Path))

        #I am obtaining the duplicated genes from each Clade (all repetitions inside each Clade are already eliminated)
        Duplications_by_Clade = []
        for i in range(0,len(Tree_Path)):
                duplication_events_table = DUPLICATIONS_DF[(DUPLICATIONS_DF['Species Tree Node'] == Tree_Path[i]) & (DUPLICATIONS_DF['Support'] >= Support)] 
                Duplications_by_Clade.append(extract_duplication_events(duplication_events_table))

        #Now I am comparing each Clade with each other in order to eliminate in between clades repetitions. As it is sorted from newer duplications to older the newest duplicated genes will be kept and the older one will be removed
        Copies = []
        for i in range(0,len(Duplications_by_Clade)):
            for z in range (1,len(Duplications_by_Clade)):
                for elemento in Duplications_by_Clade[z]:
                    if elemento in Duplications_by_Clade[i]:
                        Copies.append(elemento)

        for i in range(1,len(Duplications_by_Clade)):
                for z in range (0,len(Duplications_by_Clade[i])):
                    if Copies[z] in Duplications_by_Clade[i]:
                        Duplications_by_Clade[i].remove(Copies[z])
                    else:
                        continue

        #Here we are getting the duplication events in each node. 
        Duplication_events = []
        for i in range(0,len(Tree_Path)):
            table = DUPLICATIONS_DF[(DUPLICATIONS_DF['Species Tree Node'] == Tree_Path[i]) & (DUPLICATIONS_DF['Support'] >= Support)]
            Duplication_events.append(len(table))

        #It is important to be certain about the duplication events of the reference specie, that is why we are calculating that number in a diferent way:
        Reference_Specie_DF = DUPLICATIONS_DF[(DUPLICATIONS_DF['Species Tree Node'] == Reference_Specie) & (DUPLICATIONS_DF['Support'] >= Support)] 

        Orthogroups_copies = []
        for Orthogroup in Reference_Specie_DF['Orthogroup']:
            Orthogroups_copies.append(Orthogroup)

        Orthogroups = []    
        for Orthogroup in Orthogroups_copies:
            if Orthogroup not in Orthogroups:
                Orthogroups.append(Orthogroup)

        Total_genes_Orthogroup = []
        for Orthogroup in Orthogroups:
            Genes_Orthogroup = []
            ORTHOGROUP_DF = Reference_Specie_DF[(Reference_Specie_DF['Orthogroup'] == Orthogroup)]
            for i in range(0,len(ORTHOGROUP_DF)):
                Genes_Orthogroup += ((ORTHOGROUP_DF.iloc[i]['Genes 1'].replace(' ','') + ',' + ORTHOGROUP_DF.iloc[i]['Genes 2'].replace(' ','')).split(','))
            Total_genes_Orthogroup.append(Genes_Orthogroup) 

        Non_duplicated_genes_Orthogroup = []
        for element in Total_genes_Orthogroup:
            Non_duplicated_genes_Orthogroup.append(list(set(element)))

        Duplication_Events_Ref_Spe = 0
        for element in Non_duplicated_genes_Orthogroup:
            Duplication_Events_Ref_Spe += (len(element)-1)        

        #In order to normalize the duplication events we need to now the branch length of each node.
        Branch_Length = {}
        for node in S_tree.find_clades(branch_length = True):
            Branch_Length[node.name] = node.branch_length

###################################################################################################

        #Creating a tsv file to save all the duplications of each node
        Clade_Duplications = {}
        for i in range(0,len(Tree_Path)):
            if '_proteins' in Tree_Path[i]:
                Clade_Duplications[Tree_Path[i].replace('_proteins', '')] = Duplications_by_Clade[i]
            else:
                Clade_Duplications[Tree_Path[i]] = Duplications_by_Clade[i]

        Duplications_file = open(str(date.today()) + '-' + Output + "_Duplication_Analysis.tsv", "w")

        writer = csv.writer(Duplications_file, delimiter = '\t')
        writer.writerow(['Node','Genes'])
        for key, value in Clade_Duplications.items():
            for gene in value:
                writer.writerow([key, gene])

        Duplications_file.close()

###################################################################################################

        #Creating a tsv file to save the overall information of the Gene Duplications Analysis
        Duplications_Overall = open(str(date.today()) + '-' + Output + "_Duplication_Overall.tsv", "w")

        writer = csv.writer(Duplications_Overall, delimiter = '\t')

        counter = 0
        writer.writerow(['Node','Specie','Duplications', 'Dupl Events', 'Branch length', 'Normalized events'])
        for key, value in Clade_Duplications.items():
            if Reference_Specie_short in key:
                writer.writerow([key, list(Clade_Path.keys())[counter], len(value),Duplication_Events_Ref_Spe, Branch_Length.get(key + '_proteins'),(Duplication_Events_Ref_Spe)/(Branch_Length.get(key+'_proteins')*1000)])  
            if key not in Branch_Length.keys() and Reference_Specie_short not in key:
                continue
            elif Reference_Specie_short not in key:
                writer.writerow([key, list(Clade_Path.keys())[counter],len(value),Duplication_events[counter], Branch_Length.get(key),(Duplication_events[counter])/(Branch_Length.get(key)*1000)])

            counter += 1

        Duplications_Overall.close()

        print("Analysis done, new file: " + str(date.today()) + '-' + Output + "_Duplication_Analysis.tsv")
        print("Analysis done, new file: " + str(date.today()) + '-' + Output + "_Duplication_Overall.tsv")
        
        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), ' : Done duplicated genes')
        
        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), ' : Running De novo Analysis')
        
###################################################################################################        
#------------------------------------------------------------------------------------------------------------#
#DE NOVO
#------------------------------------------------------------------------------------------------------------#
        #First we transform all the files in the Orthologues folder into DataFrames
        Orthologues_Files = glob.glob(Orthologues_Folder)

        Orthologues_Tables = []
        for file in Orthologues_Files:
            data = pd.read_csv(file, sep = '\t')
            Orthologues_Tables.append(data)

        #Here we are obtaining all the orthologues in the dataframe and saving them in lists.       
        Orthologues = []
        for table in Orthologues_Tables:
            Group_genes=[]
            for gene in table[Reference_Specie]:
                Group_genes.append(gene.replace(' ',''))
            Total_genes = ','.join(Group_genes)
            Group_genes = Total_genes.split(',')
            Orthologues.append(Group_genes)

        #Species names from the Orthologues_Folder
        Species_Names = []
        for table in Orthologues_Tables:
            Species_Names.append(table.columns[2])
        Species_Names=short_names(Species_Names)

        #Para los siguientes pasos necesitas que está lista no contenga la especie referencia
        Tree_Path.remove(Reference_Specie)

        #Here we are creating a dictionary in order to now wich Orthologues are from what specie. We are interested in getting the lists of Orthologues ordered according to the tree path
        Species_Genes = {}
        for i in range(0,len(Species_Names)):
            Species_Genes[Species_Names[i]] = Orthologues[i]

        Orthologues_ordered_list = []
        for i in range(0,len(Clades_Only)):
            Orthologues_ordered_list += Species_Genes.get(Clades_Only[i])

        #Once we have all the orthologues we are interested in knowing wich ones are exclusive of the reference specie. First we have to get all the proteins used in the study.
        All_genes=[]
        for record in SeqIO.parse(Proteins_Fasta, "fasta"):
            All_genes.append(record.id)

        #Here we compare all the proteins used with those we have inside an Orthogroup and as a result we have the proteins exclusive of the reference specie.
        ReferenceSpecie_de_novo = []
        for genes in All_genes:
            if genes not in Orthologues_ordered_list:
                ReferenceSpecie_de_novo.append(genes)

        #Once the reference specie is done we have to sort de novo genes by specie.   
        nonExclusive_de_novo = []
        counter = 0
        for i in range(0,len(Clades_Only)):
            Species_de_novo = []
            counter += len(Species_Genes.get(Clades_Only[i]))
            for gene in Species_Genes.get(Clades_Only[i]):
                if gene not in Orthologues_ordered_list[counter:]:
                    Species_de_novo.append(gene)
            nonExclusive_de_novo.append(Species_de_novo) 

        De_novo_genes = {}
        De_novo_genes[Reference_Specie_short] = ReferenceSpecie_de_novo
        for i in range(0,len(Clades_Only)):
            De_novo_genes[Clades_Only[i]] = nonExclusive_de_novo[i]

        #Now we want to have every de novo protein idexed by specie and by Orthogroup. In order to do that we need all the Orthogroups provided by OrthoFinder
        Orthogroups_non_exclusive = pd.read_csv(Orthogroups_File, delimiter = "\t", header=0)
        Orthogroups_Exclusive = pd.read_csv(Orthogroup_UnassignedGenes, delimiter = "\t", header=0)

        #Here we are using orthogroups_todict to get each protein related to an Orthogroup.
        Orthogroups_non_exclusive_dict = orthogroups_todict(Orthogroups_non_exclusive, Reference_Specie)
        Orthogroups_exclusive_dict = orthogroups_todict(Orthogroups_Exclusive, Reference_Specie)

        #In order to have all the information of the protein in just one object we are going to create a dictionary (key = protein : values = [Orthogroup, Specie, Node])
        De_novo_gene_information = {}
        for specie in De_novo_genes:
            for gene in De_novo_genes.get(specie):
                if gene in list(Orthogroups_non_exclusive_dict.keys()):
                    De_novo_gene_information[gene] = [Orthogroups_non_exclusive_dict.get(gene), specie]
                elif gene in list(Orthogroups_exclusive_dict.keys()):
                    De_novo_gene_information[gene] = [Orthogroups_exclusive_dict.get(gene), specie]

        for element in Clades:
            if len(element) > 1:
                for specie in De_novo_gene_information:
                    if De_novo_gene_information[specie][1] in short_names(element):
                        De_novo_gene_information[specie][1] = element[1].replace('_proteins', '') + '_group'

        for element in De_novo_gene_information:
            De_novo_gene_information[element].append(Clade_Path[De_novo_gene_information[element][1]])

        #In order to obtain the number of 'de novo' events we are going to count the number of Orthogroups    
        Species_Path = list(reversed(Species_Path))
        Species_Path.append(Reference_Specie_short)
        Species_Path = list(reversed(Species_Path))

        count_Orthogroups = []
        for i in range(0, len(Species_Path)):
            Orthogroup = []
            for gene in De_novo_gene_information:
                if De_novo_gene_information[gene][1] == Species_Path[i]:
                    if De_novo_gene_information[gene][0] not in Orthogroup:
                        Orthogroup.append(De_novo_gene_information[gene][0])
            count_Orthogroups.append(len(Orthogroup))

###################################################################################################           
        #As an output we are going to get two files: one showing all "de novo" proteins and their information and another showing the sum up of this information. 
        #Creating a tsv file to save all the duplications of each node
        De_novo_file = open(str(date.today()) + '-' + Output + "_De_novo_Analysis.tsv", "w")

        writer = csv.writer(De_novo_file, delimiter = '\t')
        writer.writerow(['Genes','Orthogroup', 'Clade', 'Node'])
        for key, value in De_novo_gene_information.items():
            writer.writerow([key, value[0], value[1], value[2]])

        De_novo_file.close()
        
###################################################################################################

        #File containing the overall 'de novo' information
        De_novo_Overall = open(str(date.today()) + '-' + Output + "_De_novo_Overall.tsv", "w")

        writer = csv.writer(De_novo_Overall, delimiter = '\t')
        writer.writerow(['Node','Specie','De novo genes', 'De novo events', 'Branch length', 'Normalized events'])

        i = 0
        for specie in Clade_Path:
            #El counter te da el número de genes de novo
            counter = 0 
            for element in De_novo_gene_information:
                if De_novo_gene_information[element][1] == specie:
                    counter += 1  
            if Reference_Specie_short in Species_Path[i]:
                writer.writerow([specie,specie,counter,count_Orthogroups[i],Branch_Length.get(Clade_Path.get(Species_Path[i])+'_proteins'), count_Orthogroups[i]/(Branch_Length.get(Clade_Path.get(Species_Path[i])+'_proteins')*1000)])
            if Branch_Length.get(Clade_Path.get(Species_Path[i])) is None and i != 0:
                continue
            elif Branch_Length.get(Clade_Path.get(Species_Path[i])) is not None:
                writer.writerow([Tree_Path[i-1],specie,counter,count_Orthogroups[i],Branch_Length.get(Clade_Path.get(Species_Path[i])),count_Orthogroups[i]/(Branch_Length.get(Clade_Path.get(Species_Path[i]))*1000)])
            i += 1
        De_novo_Overall.close()    

        print('Analysis done, new file: ', str(date.today()) + '-' + Output + "_De_novo_Analysis.tsv")
        print('Analysis done, new file: ', str(date.today()) + '-' + Output + "_De_novo_Overall.tsv")
        
###################################################################################################

        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), ' : Done de novo genes')

        print('Analysis finished')

    else:
        Configurations = ['Reference_Specie:', 'Duplications.tsv:', 'Species_Tree:', 'Non_labeled_Species_tree:', 'Orthologues_folder:', 'Proteins_fasta:', 'Orthogroups_file:', 'Orthogroup_UnassignedGenes.tsv:', 'Output:', 'Support:']
        with open ('settings.txt','w') as file:
            file.writelines("%s\n" % i for i in Configurations)
        print('Settings file was not found')
        print('settings.txt file created, needs to be filled')