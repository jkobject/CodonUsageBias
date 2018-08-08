# METADATA INFO
- the molecular property of each amino acid molecules
- the genomes of Yun's dataset with their links
- the name of the genomes ordered by their ancestry line
- the name of the genomes with their ftp server links to ensembl db
- an example of GFF3 files that can be downloaded for each species
- the file metada for the kingdom's available species
- the mean abundance in proteins per cells for cerevisiae

## protdata/tob_currated.csv
Entries in the curated dataset are colour-coded with reference to the steps shown below.

with Gene Name	Molecular Weight (Da)	Length (Amino Acids)	Protein Abundance (molecules per cell)	Protein decay rate (min-1)	mRNA Abundance (molecules per cell)

### Protein abundances:
Studies were scaled to a common mean as described in the text. The average of all available scaled values was then used to create the curated dataset. (4242 entries).
All Proteins which could be successfully tagged with both GFP and TAP-tags, but where expression could be observed with neither tag, and where no published evidence to the contrary could be located, where assigned an abundance value of 0 (1333 entries).
For all remaining genes for which mRNA abundance values were available, the average ratio of experimental protein to mRNA abundances (3450 proteins/mRNA) was used to predict protein abundance (450 entries).
551 genes remain without assigned protein abundance.

### mRNA abundances:
Five different transcriptome studies were scaled to match their mean mRNA abundance, and averages were calculated from values in the scaled studies (6025 entries).
551 entries remain without assigned mRNA abundance.

### Protein half-lives:
Where available, data from the genome wide study by Belle et al. were used in the curated dataset. Where no data were available from Belle et al. but data were reported in other studies, these latter values were used (3178 entries).
For proteins with reported half-lives of >300 minutes, or where the literature only reported "no decay", the decay rate was set to 0 (388 entries).
For all other entries, the average decay rate of all experimentally determined entries was used as "best guess" (2760 entries).

## protdata/Amino Acid Properties.csv
Data sources for physical properties of amino acids:

### General sources:

Molecular_Weight - the molecular weight of the amino acid without crystal water

NH2_pKA,COOH_pKA,sidechain_pKA - the dissociation constant for acid/ base groups

pI - the isoelectric point

no_atoms - the number of atoms making up the amino acid

>>>Zamyatnin, A.A. (1972) Prog. Biophys. Mol. Biol., 24: 107-123

volume - solution volume in AA^3

>>>Monera et al., J. Protein Sci. 1: 319-329 (1995).

Hydrophobicity_index - a relative measure of the hydrophobicity 

>>>This study (with data from Henikoff and Henikoff (1992) Proc. Natl. Acad. Sci. USA 89: 10915–10919. 

Conservation_Index: calulated from the BLOSUM62 substitution matrix as abs(sum(AA(non-self))/AA(self))

>>>Barton et al. (2010) PLoS One 5: e11935

rel_C_cost, rel_N_cost, rel_S_cost; cost in terms of individual elements
rel_glucose; relative synthesis cost in terms of glucose equivalents (A_glucose in paper)

>>>Craig and Weber (1998) Molecular Biology and Evolution 15: 774–776.

synthesis_steps, number of steps required between central metabolism and fully formed amino acid

##ENSEMBLMETA

JSON format file containing full details of the analyses available for each genome as an array named "genome". Each element of the array contains the following keys:
 - name - full name of species/strain
 - species - computationally safe version of species as used internally by Ensembl
 - division - division of Ensembl Genomes (e.g. EnsemblFungi)
 - taxonomy_id - NCBI taxonomy ID
 - assembly_name - version of assembly
 - genebuild_id - version of genebuild
 - pan_species - 0/1 indicating if the genome is included in the pan taxonomic compara
 - db_name - name of MySQL core database
 - species_id - ID of genome with MySQL core database
 - annotation - associated array containing the following keys:
 - nProteinCoding - number of protein-coding genes in this genome
 - nInterPro - number of InterPro entries annotated to protein-coding genes in this genome
 - nInterProDomains - number of InterPro domains annotated to protein-coding genes in this genome
 - nGO - number of GO terms annotated to protein-coding genes in this genome
 - nUniProtKBSwissProt - number of UniProtKB/Swiss-Prot entries annotated to protein-coding genes in this genome
 - nUniProtKBTrEMBL - number of UniProtKB/TrEMBL entries annotated to protein-coding genes in this genome
 - nProteinCodingUniProtKB - number of protein-coding genes with UniProtKB cross-references
 - nProteinCodingUniProtKBSwissProt - number of protein-coding genes with UniProtKB/Swiss-Prot cross-references
 - nProteinCodingUniProtKBTrEMBL - number of protein-coding genes with UniProtKB/TrEMBL cross-references
 - nProteinCodingInterPro - number of protein-coding genes with InterPro domains
 - nProteinCodingGO - number of protein-coding genes with GO terms
 - variation - associative array containing the following keys:
 - variations - counts of variations keyed by source names
 - structural_variations - counts of structural_variations keyed by source names
 - phenotypes - names of phenotypes with counts of associated variations
 - genotypes - names of samples with counts of associated non-reference variations
 - compara - associative array where keys are names of comparative analyses that the genome is involved in, and values are lists of genomes included in those analyses
 - features - associative array with the following keys
 - proteinAlignFeatures - associative array of analysis names and counts of protein-genome alignments
 - dnaAlignFeatures - associative array of analysis names and counts of DNA-genome alignments
 - repeateFeatures - associative array of analysis names and counts of repetitive sequence features
 - simpleFeatures - associative array of analysis names and counts of other genomic features
 - bam - associative array where keys are names of BAM file types (dna, rnaseq etc.), and values are arrays of associative arrays, one per BAM file with the following keys:
 - source_name - descriptive name of BAM file
 - source_url - URL where BAM file can be found
 - id - unique identifier of BAM file
-- description - longer description of the BAM file

## protdata/PDI_substrates

## Metaphylo

## homolist.json

## short Homolist.json

## 3Dmodels

## names_with_links

## order_name461

## Yun_Species_Context