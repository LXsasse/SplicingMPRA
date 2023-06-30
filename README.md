# Impact of random 50-base sequences inserted into an intron on splicing in Saccharomyces cerevisiae 

Code for some of the data analysis performed in the [paper](https://www.biorxiv.org/content/10.1101/2023.06.21.545966v1)

Linear models were trained as defined in `run_linear.sh`

Data processing was performed as in `run_data_processing.sh`

CNN models were trained as described in `run_cnns.sh`

## Results
To determine sequence patterns in the random mRNA element (N50) that influence splicing efficiency (SE), we applied five different convolutional neural network (CNN) architectures with sequence-structure inputs of various lengths (Methods). To train  and evaluate these models, we selected 140,017 sequence-SE pairs that showed high concordance between replicates. We trained all models on the mean efficiency from the two replicates. We trained on 80% of the sequence-SE pairs, validated on another 10%, and tested performance on the left-out 10%. 

None of the tested models was able to explain more than 7% of the variance that was measured in the MPRA (i.e. max Pearson’s R = 0.265, Table 1). Downsampling of the number of sequences in the training set suggested that the CNN models would require much more training sequences to learn the patterns that are controlling splicing efficiency in the MPRA. 

We also applied linear models that used only 15 pre-selected features describing splicing, and experiment related sequence features, as well as structural descriptors of the N50. These linear models performed almost on par with the way more complex CNN models (Table 2). The most predictive feature of the linear models was the predicted minimum free energy of the structural ensemble (MFEE), followed by the presence of an additional 5’ splice site, the presence of GU-rich hexamers that were enriched in sequences with low splicing efficiency, and the GC-content of the N50. The lower the MFEE, i.e. the more stable the structure of the transcript, the lower was the splicing efficiency, suggesting that formation of strong RNA structures influences splicing efficiency, likely by controlling access of splicing factors to splice sites and branch points. This effect of the RNA structure on splicing is likely the reason why our CNN models require so much more data to achieve higher correlations. While CNNs are great at learning combinations of sequence motifs, for example from RNA binding proteins, they require more data to learn complex non-linear interactions that are represented in RNA structure, where a multitude of sequences in different arrangements can form the same secondary or tertiary structures. On the other hand, current RNA structure prediction methods have limited accuracy, and it is unclear if the structural features that we extracted fully captured the essential properties that influence splicing. 

Lastly, we trained an elastic-net model on sequence k-mers within the N50 to determine which sequence patterns have the strongest influence SE. These models predicted GUG elements to negatively affect SE while elements with AAC elements positively affected SE. We note that the 3’UTR is enriched in CAC elements that could base-pair with the GUG-elements in the N50, and therefore could explain the reduction of SE due to lower accessibility to splice-sites. However, our RNAfold predictions did not confirm the formation of base-pairs between the N50 and the 3’UTR, suggesting that either RNAfold predictions lack accuracy or that the formation of secondary structure is insufficient to explain the splicing efficiency.

[Methods and Tables](https://docs.google.com/document/d/1yEEmDrk8YUNmIng5fCt_-F0GEAKnifXQ2I-4YX7bc4k/edit?usp=sharing)
