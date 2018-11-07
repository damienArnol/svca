

data_dir='../../examples/data/IMC_example/Cy1x7/' 

output_dir='./test_svca'
mkdir $output_dir
protein_index=23
normalisation='quantile'

python run_indiv.py $data_dir $output_dir $protein_index $normalisation 
