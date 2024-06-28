# Interpreting ML kinetics

This repository contains the data, code, and trained models for the work entitled "Kinetic predictions for SN2 reactions using the BERT architecture: Comparison and interpretation". 

<img width="1109" alt="image" src="https://github.com/C-Wils/InterpretingMLKinetics/assets/88711576/ad528e27-c8ca-4b3c-85f1-096ad9361106">

## Python versions and libraries

Code for BERT model training and evaluation, and calculation of Integrated Gradients, was tested using Python version 3.6.13 and the following libraries: pandas, numpy, scipy, rdkit, regex, openpyxl, scikit-learn, transformers, torch, datasets, dotenv, ax-platform and captum (versions specified in Scripts/Training_and_evaluation/BERT/requirements.txt)

All other code was tested using Python version 3.9.12 and the following libraries: pandas, numpy, scipy, matplotlib,  rdkit, regex, scikit-learn, CIMtools, CGRtools and treeinterpreter, pywaffle (versions specified in ./requirements.txt)

