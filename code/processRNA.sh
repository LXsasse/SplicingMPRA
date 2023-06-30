# Run RNAfold to generate <seq>dp.ps files in a directory 
# download RNAfold from and follow installation instructions: https://github.com/ViennaRNA/ViennaRNA
RNAfold -i <filename>.fasta -g -p --noPS > <filename>.struc.out

codedir=/path_to_repsitory/code/

# In case one wants to plot RNA strucutre matrix
python3 ${codedir}/plot_ps_matrix.py seq_AGAACGGAACCGGGC_dp.ps

# Transform ps files into matrices and attach as features to sequence representations.
#infile=$1
#awk -v var="$i" '{if ($4=="ubox") print $1,$2,$3,var}'
python3 ${codedir}/transformpstoMatrix.py input_folder/with/psfiles
# potentially get energy parameters of the bp and multiply them by the bp-probabilities 

# Transform output fasta into pure fasta
python3 ${codedir}/output_to_fasta.py <filename>.struc.out

# Use transformed fasta that only has sequence and bracket notion of MFE structure to get secondary structure code
#### First, download forgi, install and use adopted rnaConvert.py file to convert minimum free energy structure into structure element representation. 
# https://viennarna.github.io/forgi/
python3 ${codedir}/rnaConvert.adopt.py -T element_string <filename>.struc.fasta > filename>.element_struc.fasta

# Use secondary structure code file and fasta to get structure representations
python3 ${codedir}/multistructure_representation.py N50_extended.struc.fasta N50_extended.element_struc.fasta 


