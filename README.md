# Percolator with FDR control

This contains the code for the associated manuscript, titled "Splitting the decoys into training and estimating: how to train a Percolator-like post-processor and maintain control of the false discovery rate". 

The important folder is *functions* which contains the actual scripts that implements our procedure using Python, from the command line. The main python script is the file, do_FDR_percolator.py, and the auxiliary functions can be found under percolator_functions.py and utility_functions.py. This python implementation may be used as interim until the development in C (Percolator) or the mokapot module is updated. 

A note to Lukas: for the development in C, please refer to the link [here](https://unisydneyedu-my.sharepoint.com/:f:/g/personal/jfre0619_uni_sydney_edu_au/Et239lX9DY9LjFA4GtfyHvQBMDRZbCc8dcU2QDaE_2mDPw?e=gUJDXR) which contains test files and their associated results using the current Python implementation.

