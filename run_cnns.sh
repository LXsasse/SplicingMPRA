#!/bin/bash

##### Download code from https://github.com/LXsasse/DRG/ ######

infile=$1
jobid=$2

cv=0
cvfile='N50_cv10.txt'

bs=16
lr='1e-4'
patience=40
epochs=1000

outfunc='Linear'

NOW=$( date '+%F_%H-%M-%S' )

ctrl=Models/${jobid}_out${NOW}.txt
echo $NOW > $ctrl
echo 'INFILE : '${infile} >> $ctrl

codedir=/path_to_code_directory/DRG/

# regCNN
echo '++++++++++ regression CNN +++++++++' >> $ctrl
python ${codedir}cnn_model.py $infile N50_meanefficacy.tsv --delimiter ',' --outdir Models/ --crossvalidation $cvfile $cv 0 --maketestvalset --cnn 'loss_function=MSE+validation_loss=Correlationclass+num_kernels=200+kernel_bias=False+l_kernels=15+kernel_function=ReLU+max_pooling=False+weighted_pooling=True+pooling_size=15+outclass='${outfunc}'+shift_sequence=0+random_shift=False+epochs='${epochs}'+patience='${patience}'+batchsize='${bs}'+lr='${lr}'+write_steps=1+device=cuda:0' --save_correlation_perclass --save_mse_perclass >> $ctrl

# shallow CNN
echo '++++++++++ shallow CNN +++++++++' >> $ctrl
python ${codedir}cnn_model.py $infile N50_meanefficacy.tsv --delimiter ',' --outdir Models/ --crossvalidation $cvfile $cv 0 --maketestvalset --cnn 'loss_function=MSE+validation_loss=Correlationclass+num_kernels=200+kernel_bias=False+l_kernels=15+kernel_function=ReLU+max_pooling=False+weighted_pooling=False+pooling_size=1+dilated_convolutions=4+l_dilkernels=4+dilations=[2,4,8,16]+conv_increase=1.05+dilweighted_pooling=2+nfc_layers=2+outclass='${outfunc}'+shift_sequence=0+random_shift=False+epochs='${epochs}'+patience='${patience}'+batchsize='${bs}'+lr='${lr}'+write_steps=1+device=cuda:0+conv_batch_norm=True+fc_dropout=0.1' --save_correlation_perclass --save_mse_perclass --save_predictions >> $ctrl

# Deep CNN
echo '++++++++++ Deep CNN +++++++++' >> $ctrl
python ${codedir}cnn_model.py $infile N50_meanefficacy.tsv --delimiter ',' --outdir Models/ --crossvalidation $cvfile $cv 0 --maketestvalset --cnn 'loss_function=MSE+validation_loss=Correlationclass+num_kernels=200+kernel_bias=False+l_kernels=15+kernel_function=ReLU+max_pooling=False+weighted_pooling=False+pooling_size=1+dilated_convolutions=5+l_dilkernels=5+conv_increase=1.05+transformer_convolutions=4+trdilations=[2,4,8,16]+l_trkernels=3+trweighted_pooling=2+nfc_layers=2+outclass='${outfunc}'+shift_sequence=0+random_shift=False+epochs='${epochs}'+patience='${patience}'+batchsize='${bs}'+lr='${lr}'+write_steps=1+device=cuda:0+conv_batch_norm=True+fc_dropout=0.1' --save_correlation_perclass --save_mse_perclass >> $ctrl

alr='1e-7'
echo '++++++++++ Deep attention CNN +++++++++' >> $ctrl
python ${codedir}cnn_model.py $infile N50_meanefficacy.tsv --delimiter ',' --outdir Models/ --crossvalidation $cvfile $cv 0 --maketestvalset --cnn 'loss_function=MSE+validation_loss=Correlationclass+num_kernels=200+kernel_bias=False+l_kernels=15+kernel_function=ReLU+max_pooling=False+weighted_pooling=False+pooling_size=1+dilated_convolutions=5+l_dilkernels=5+dilations=1+conv_increase=1.05+n_attention=4+n_distattention=6+dim_embattention=196+attentionweighted_pooling=2+outclass='${outfunc}'+shift_sequence=0+random_shift=False+epochs='${epochs}'+patience='${patience}'+batchsize='${bs}'+lr='${alr}'+write_steps=1+device=cuda:0+conv_batch_norm=True+attention_dropout=0.1+fc_dropout=0.1' --save_correlation_perclass --save_mse_perclass >> $ctrl

alr='2e-8'
echo '++++++++++ Deep attention multi CNN +++++++++' >> $ctrl
python ${codedir}cnn_model.py $infile N50_meanefficacy.tsv --delimiter ',' --outdir Models/ --crossvalidation $cvfile $cv 0 --maketestvalset --cnn 'loss_function=MSE+validation_loss=Correlationclass+num_kernels=200+kernel_bias=False+l_kernels=15+kernel_function=ReLU+max_pooling=False+weighted_pooling=False+pooling_size=1+dilated_convolutions=5+l_dilkernels=5+dilations=1+conv_increase=1.05+n_attention=4+n_distattention=6+dim_embattention=196+transformer_convolutions=4+trdilations=1+l_trkernels=3+trweighted_pooling=2+outclass='${outfunc}'+shift_sequence=0+random_shift=False+epochs='${epochs}'+patience='${patience}'+batchsize='${bs}'+lr='${alr}'+write_steps=1+device=cuda:0+conv_batch_norm=True+attention_dropout=0.1+fc_dropout=0.1' --save_correlation_perclass --save_mse_perclass >> $ctrl


