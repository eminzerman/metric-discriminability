# metric-discriminability

This project consists of the Matlab code of the objective quality metric discriminability algorithm. This code is developed to understand the dicrimination power of different objective quality metrics on a subjectively-annotated quality database. It can be used as an alternative to the commonly used statistical analysis methods --which are described in ITU-T Rec. P.1401.

## Use

You can simply clone or download the project and run **metricDiscriminability** within Matlab. It requires the objective quality results, subjective scores as MOS values and Std values, number of stimuli, and number of subjects. 

The resulting Matlab struct includes the area under curve (AUC) of the receiver operating characteristic (ROC) curve. This struct also include the optimal threshold value for metric discrimination and balanced accuracy for the cases of false positive rate (FPR) < 0.05, FPR < 0.25, and the maximum accuracy case.

## Licence

This software is provided under GNU General Public Licence.

If you use this software for research purposes, we kindly ask you to cite our paper below:
* E. Zerman, G. Valenzise, and F. Dufaux. **_"An extensive performance evaluation of full-reference HDR image quality metrics."_** Quality and User Experience, volume 2, April 2017.