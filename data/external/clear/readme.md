CLEAR uses hg18 (NCBI36) in their genomic coordinates but we are using GRCh38 in our analysis. 

This data "clear_grch38_coordinates.bed" is the transformed version of CLEAR coordinates.

Methods:

- import CLEAR data to a pandas dataframe
- select only relevant (chr, start, end) positions
- use pybed submodule (found in fuc package) to transform df into a .bed file
- use ENSEMBL assembly converter to convert into GRCh38 coords.

(https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter)