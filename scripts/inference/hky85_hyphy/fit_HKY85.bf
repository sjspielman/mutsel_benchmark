LIKELIHOOD_FUNCTION_OUTPUT = 1;
DataSet ds = ReadDataFile ("hyin.txt");
DataSetFilter filtered_nuc_data = CreateFilter(ds,1);
global k;
global time;
global aux_1 = 0.5; aux_1 :< 1;
global aux_2 = 0.5; aux_2 :< 1;
global aux_3 = 0.5; aux_3 :< 1;

ML_Freqs = {{
aux_1, /*A*/
(1-aux_1)*aux_2, /*C*/
(1-aux_1)*(1-aux_2)*aux_3, /*G*/
(1-aux_1)*(1-aux_2)*(1-aux_3) /*T; note that, in this manner, the 4 frequencies will always be in [0,1] and sum to 1*/
}};
HKY85 = {{,time,k*time,time}
         {time,,time,k*time}
         {k*time,time,,time}
         {time,k*time,time,}};
         
Model HKYmodel = (HKY85, ML_Freqs, 1);
UseModel (USE_NO_MODEL);
UseModel(HKYmodel);
Tree nuc_tree = DATAFILE_TREE;
LikelihoodFunction LikFn_nuc = (filtered_nuc_data, nuc_tree);
Optimize (paramValues_nuc, LikFn_nuc);
fprintf(stdout, LikFn_nuc);