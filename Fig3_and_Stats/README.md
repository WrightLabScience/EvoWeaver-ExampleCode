# Generating Statistics 

`GenComplexStatistics.R` and `GenModuleStatistics.R` show how the files in Fig2_FigS1 are created. 
These take the component scores and generate ROC curve values, as well as train/test ensemble classifiers.

# Multiclass classification

Multiclass classification for Figure 3 is done with `CVMulticlassClassif.R`. 
The highest confidence misclassifications are shown in `HighestMisclassifications.xlsx`. 
These data are the source for Figure 3 and Figure S2, though Biorender was used to improve the visuals of 
the raw outputs. 

# Note
Some datafiles are too large to include on GitHub; see Zenodo for these files.
