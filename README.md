# Time-series-transformers (TST) for Crystallization Systems
Novel time-series-transformer framework that can show model predictions across N different chemical systems (e.g., crystallization), providing high system-to-system (S2S) interchangeability and more accurate time-series predictions than current SOTA ML models (e.g., RNN, LSTM, CNN, etc.). 


Project description:
 
This project aimed to develop a novel TST model that can show model predictions across 20 different (or N different) crystal systems and be better than current SOTA ML models (e.g., RNN, LSTM, CNN, etc.) 
 
1. The initial crystallization simulation code is included in ‘Crystallization_case_study.PDF’ file and ‘CASE_1_Dextrose.ipynb’ file.
 
2. Parallelization on high-performance computing (HPC) was performed for “N” different crystallization systems, and this is included in the folder ‘Generalized_crys_simulations’, which generates data of ~15M+ datapoints, and also compiles them in an appropriate .H5 (memory-efficient) file. It has many different scripts that need to be included for properly running as well. 
 
3. The basics and explanation of a TST framework is included in the ‘TST_base_code.PDF’ file and ‘TST_base_code.ipynb’ file.
 
4. After training model testing results and validation are shown in ‘TST_model_testing.PDF’. This requires the actual model (a PyTorch model to be installed), and all the various dataset files. Thus, python code for this part is not included. But the .PDF file should give all the necessary information as that is what I use for testing purposes.
 
5. To get an overall summary of the work, I have attached a recent paper we published on this work, and I have included a .PPT file to get an overview of the current work (which is still under review in another journal). 

