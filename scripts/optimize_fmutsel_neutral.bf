_Genetic_Code = {
                     {14,13,14,13,7,7,7,7,19, 5,19, 5,2,2,3,2,
                      12,11,12,11,6,6,6,6,19,19,19,19,1,1,1,1,
                      16,15,16,15,8,8,8,8,20,20,20,20,4,4,4,4,
                      10,9, 10,9, 5,5,5,5,10,17,18,17,1,0,1,0}
                };
function buildCodonFrequencies1x4(obsF)
{
    PIStop = 1.0;
    result = {61, 1};
    hshift = 0;

    for (h=0; h<64; h=h+1)
    {
        first = h$16;
        second = h%16$4;
        third = h%4;
        if (_Genetic_Code[h]==10) 
        {
            hshift = hshift+1;
            PIStop = PIStop - obsF[first]*obsF[second]*obsF[third];
            continue; 
        }
        result[h-hshift][0]=obsF[first]*obsF[second]*obsF[third];                
    }
    return result*(1.0/(PIStop));
}


OPTIMIZATION_PRECISION = 0.000001;
DataSet ds = ReadDataFile ("hyin.txt");
DataSetFilter filtered_nuc_data = CreateFilter(ds,1);
global k;
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
//fprintf(stdout, paramValues_nuc);

//fprintf(stdout, "\nNucleotide model optimized. Now codon model\n.")

/* Now that we have frequencies and kappa, optimize branch lengths with FMutSel neutral model (all substitutions equally probable, so all selection coefficients == 1 and evolution dictated by mutation alone). */


global freqA := aux_1;
global freqC := (1-aux_1)*aux_2;
global freqG := (1-aux_1)*(1-aux_2)*aux_3;
global freqT := (1-aux_1)*(1-aux_2)*(1-aux_3);
#include "fmutsel_neutral.mdl";


f1x4 = buildCodonFrequencies1x4({{freqA, freqC, freqG, freqT}}); // F1x4
DataSet ds2 = ReadDataFile ("hyin.txt");
DataSetFilter filtered_codon_data = CreateFilter(ds2, 3, "", "", "TAA,TAG,TGA");

Model FMutsel_Neutral = (fmutsel_neutral, f1x4, 0);
UseModel (USE_NO_MODEL);
UseModel(FMutsel_Neutral);
Tree codon_tree = DATAFILE_TREE;
ReplicateConstraint("this1.?.k:=this2.?.k__", nuc_tree, codon_tree); // save the kappa from before!
LikelihoodFunction LikFn_codon = (filtered_codon_data, codon_tree);
Optimize (paramValues_codon, LikFn_codon);
fprintf(stdout, LikFn_codon);


