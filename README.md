# FitnessLandscape_SD

<p align=center>
<a target="_blank"><img src="https://img.shields.io/powershellgallery/p/Az.Storage.svg"></a>
<a target="_blank" href="https://www.python.org/downloads/" title="Python version"><img src="https://img.shields.io/badge/python-%3E=_3.6-green.svg"></a>
<a target="_blank" href="https://www.python.org/downloads/" title="Python version"><img src="https://img.shields.io/badge/Julia-v1.1-blueviolet.svg"></a>
<a target="_blank" href="https://opensource.org/licenses/MIT" title="License: MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg"></a>
</p>


This is the repository to store the  code for implementing the analysis and plotting in publication "Global fitness landscapes of the Shine-Dalgarno sequence".

# Dependencies
* Python 3.6
  * packages required: pandas, matplotlib, numpy, itertools, scipy,sklearn
* R
  * packages required: lattice, vioplot
* Julia Language (>= v1.0)
  * packages required: DelimitedFiles, Plots, Random, Statistics, Distributions
* Matlab 

# Usage
To run the codes and check the results, please see the next section to find the designated main file.
Take python as example: 
'''
git clone https://github.com/skyshine102/FitnessLandscape_SD.git
python ./FitnessLandscape_SD/src/Python/main_all_script.py
'''  
Note that you need to configurate the path within the code.

# Table of Contents

### Python 
  * Fitness Distribution
  * Context dependence of the G-P associations
  * Relationship between the finess of a variant and the mean fitness of its single-mutation neighbors
  * Relationship between the fitness of a SD variants and the abundance of beneficial, neutral, and deleterious mutations
  * Relationship between fitness and nucleotide composition
  * Relationship between the nucleotide content and fitness
  * Effect of single nucleotides on fitness
  * Effect of pairwise epistasis on fitness
  * Explanatory power (R^2) of single nucleotide effects and pairwise epistasis

### Julia 
  * Relationship between fitness and sequence identity

### R
  * Mutation Effect
    * Dependence of G &rarr; C mutational effects on nucleotide positions and the fitness of genetic backgrounds
    * mean trend + detailed plots
  * Correlation between the SD:aSD duplex length, base-pairing energy and fitness

# Contact Information
For more info. about the project, please contact
* PI: David H. H. Chou (chouhh@ntu.edu.tw)

As for technical support, feel free to contact
* Jeremy Jahn (skyshine102@gmail.com)
* Antony Kuo (b04b01066@ntu.edu.tw)

