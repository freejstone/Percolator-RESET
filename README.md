# Percolator RESET

Percolator RESET (REScoring via Estimating and Training) is a post-processor used to increase the number of peptide discoveries in MS/MS database searches. Percolator RESET implements the Percolator algorithm while maintaining rigorous theoretical false discovery rate (FDR) control. They key idea is that, rather than using the entire set of decoys wins for both training and estimating, a random subset of decoys is chosen for each task. One subset is used to train the Percolator algorithm which utilizes a host of feature information to separate correctly identified peptides from incorrect ones, while the other decoy set is used for controlling the FDR.

This repository contains the code for the associated manuscript, titled "[How to train a post-processor for tandem mass spectrometry proteomics database search while maintaining control of the false discovery rate](https://www.biorxiv.org/content/10.1101/2023.10.26.564068v1)". 

The important folder is *percolator_RESET* which contains the actual scripts that implements our procedure using Python, from the command line. This python implementation may be used as interim until the development in C (Percolator). For further details on implementation, see the wiki page.

A note to Lukas: for the development in C, please refer to the link [here](https://unisydneyedu-my.sharepoint.com/:f:/g/personal/jfre0619_uni_sydney_edu_au/Et239lX9DY9LjFA4GtfyHvQBMDRZbCc8dcU2QDaE_2mDPw?e=gUJDXR) which contains test files and their associated results using the current Python implementation.

For installation, download the latest release from Github, and install using pip in the same directory as setup.py using `pip install .`
