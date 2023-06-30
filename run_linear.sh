# Run feature based models to determine feature importance and as a baseline for more complex CNNs.
inputs='Seqfeat_input.npz Seqfeatkmer3_input.npz' 
kmers='Seqkmer3_input.npz Seqkmer4_input.npz Seqkmer5_input.npz Seqkmer6_input.npz'
outputs='Pred_meanefficiency.txt'
cvfile='N50_cv10.txt'
for ip in $inputs
do
	for i in {0..9}
python FeaturePredictor.py $ip $outputs --crossvalidation --crossvalidation $cvfile $i --combine_train_and_val --model_params LR 
python FeaturePredictor.py $ip $outputs --crossvalidation --crossvalidation $cvfile $i --combine_train_and_val --model_params RF
done

for k in $kmers
do
python FeaturePredictor.py $k $outputs --crossvalidation $cvfile 0 --combine_train_and_val --model_params EN alpha=${al}+l1_ratio=0.9+warm_start=True+max_iter=5000
done



