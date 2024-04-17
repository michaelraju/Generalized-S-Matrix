### About this repository 
Code packages [![DOI](https://zenodo.org/badge/785200391.svg)](https://zenodo.org/doi/10.5281/zenodo.10963077)

Saved dataset [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10960970.svg)](https://doi.org/10.5281/zenodo.10960970)

Manuscript preprint [![arXiv](https://img.shields.io/badge/arXiv-2404.10119-b31b1b.svg)](http://arxiv.org/abs/2404.10119)

This repository hosts the Matlab based code packages associated with our manuscript titled “Modeling scattering matrix containing evanescent modes for 
wavefront shaping applications in disordered media”, Michael Raju, Baptiste Jayet and Stefan Andersson-Engels, Tyndall National Institute, University College Cork, Ireland.
 In addition to the main paper, there is a supplementary 
document and a user manual. The [user manual](/User_manual.pdf)
 is provided here within the repository. It is recommended to go through the user manual first to understand the directory 
structure of the code packages. As explained in the user manual, there is a Zenodo repository with [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10960970.svg)](https://doi.org/10.5281/zenodo.10960970)
 which hosts the saved run-data. One may download the saved
run data (saved-data-Code-Package-1.mat and saved-data-Code-Package-2.mat) from 
the Zenodo repository and add it to the code package folders (Code-Package-1 and Code-Package-2) respectively. Loading the saved run data
(by setting *new\_run\_flag=0* in the main.m file) 
helps to visualize the results presented in the paper, without actually performing a new computational run from scratch. On the other hand, setting the flag 
*new\_run\_flag=1* yields a fresh computational run, initializing a new disorder.  



### Contribution
The code packages and the associated analytical and numerical formulations were developed by Michael Raju as part of his [PhD thesis](https://hdl.handle.net/10468/14107).
Baptiste Jayet and Prof. Stefan Andersson-Engels were involved in the PhD supervision. 

### Acknowledgments
We would like to thank Science Foundation Ireland (SFI) for funding the research through Stefan's professorship grant *``Novel applications and techniques for in-vivo
optical imaging and spectroscopy”* (SFI/15/RP/2828 and SFI/22/RP-2TF/10293).
