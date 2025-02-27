<img src="./figs/logo.png" alt="coremof2024" width="500"> 

#### :sparkles: Download CoRE MOF Database in [here](https://zenodo.org/uploads/14216942)                                                             
                                                                                                                 
This repository includes tools developed to collect, curate, and classify Computation-Ready, Experimental MOF database. 

**Developed by:** [Guobin Zhao](https://github.com/sxm13)                                                                                                                                  

[Web Application](https://mof-db.pusan.ac.kr/) for Database and Tutorial.                     

### Direcotry organization                                
- 0-database: This folder includes script to check publication dates and publisher infomration from a DOI input.
- 1-clean:
  -   preprocess.py: a pre-processing script to split mutil CIFs and making them primitive and P1. 
  -   clean.py: a solvent removal script to remove free and bound solvents <sup>[1](),[2](https://doi.org/10.1021/acs.jced.9b00835),[3](https://doi.org/10.1021/cm502594j)</sup>. 
- 2-NCRCheck:
  -   check MOFs by Chen and Manz<sup>[1](https://doi.org/10.1039/D0RA02498H),[2](https://doi.org/10.1039/C9RA07327B)</sup>
  -   [mofchecker](https://github.com/kjappelbaum/mofchecker)
- 3-MOFFeatures:
  - pore properties by [Zeo++](https://github.com/richardjgowers/zeoplusplus)<sup>[1](https://doi.org/10.1016/j.micromeso.2011.08.020),[2](https://doi.org/10.1021/ci200386x)</sup>, [open metal site](https://github.com/emmhald/open_metal_detector)<sup>[1](https://doi.org/10.1021/acs.jced.9b00835)<sup>
  - topology by [CrystalNets](https://github.com/coudertlab/CrystalNets.jl)<sup>[1](https://doi.org/10.21468/SciPostChem.1.2.005)<sup>
  - [revised autocorrelation functions](https://github.com/hjkgrp/molSimplify)<sup>[1](https://doi.org/10.1002/jcc.24437),[2](https://doi.org/10.1021/acs.iecr.8b04015)<sup>
  - others (mass, space group)          
- 4-MLFeatures:
  - predict partial atmoic charges by [PACMAN](https://github.com/mtap-research/PACMAN-charge)<sup>[1](https://doi.org/10.1021/acs.jctc.4c00434)<sup>
  - [heat capacity](https://github.com/SeyedMohamadMoosavi/tools-cp-porousmat)<sup>[1](https://doi.org/10.1038/s41563-022-01374-3)<sup>
  - [water stability](https://zenodo.org/records/12110918)<sup>[1](https://doi.org/10.1021/jacs.4c05879)<sup>
  - [thermal and activation](https://pubs.acs.org/doi/suppl/10.1021/jacs.1c07217/suppl_file/ja1c07217_si_002.zip)<sup>[1](https://doi.org/10.1021/jacs.1c07217)<sup>
- 5-raspaTools: generate / submit input files, check jobs and get output files for
  - [RASPA2](https://github.com/iRASPA/RASPA2)<sup>[1](https://doi.org/10.1080/08927022.2015.1010082),[2](https://doi.org/10.1080/08927022.2013.819102),[3](https://doi.org/10.1080/08927022.2013.819102),[4](https://doi.org/10.1002/adts.201900135)
  - [gRASPA](https://github.com/snurr-group/gRASPA)<sup>[1](https://doi.org/10.1021/acs.jctc.4c01058)<sup>                                                                                      
- 6-others:
  - match structures
  - DDEC files to CIF
  - check that the ASR and FSR are exactly the same structure

- 7-data4API:
  - name list of CSD CIFs, ZIP for CIFs from SI
  - information for each CR / NCR structure
  - jupyter notebook examples of usage of CoRE MOF API
    
If you encounter any issues related to the CoRE MOF Database or the scripts mentioned above, please create an **[Issues](https://github.com/mtap-research/CoRE-MOF-Tools/issues/new/choose)** and submit it. We will address it as soon as possible.
