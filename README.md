# DDGun


#### INTRODUCTION
  
  DDGun is an untrained method for predicting the stability change upon mutation


#### LICENSE

  Copyright (C) 2019  Ludovica Montanucci, Emidio Capriotti 
                      and Piero Fariselli

  This program and all program in this package are free software;
  you can redistribute it and/or modify it under the terms of the
  GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.


#### CITATION

  *Montanucci L, Capriotti E, Frank Y, Ben-Tal N, Fariselli P.* (2019).
  DDGun: an untrained method for the prediction of protein stability
  changes upon single and multiple point variations.
  **BMC Bioinformatics**. 20 (Suppl 14): 335. [PMID:31266447](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6606456/pdf/12859_2019_Article_2923.pdf).  
  Supplementary Files are available at:
  [https://doi.org/10.5281/zenodo.4613881](https://doi.org/10.5281/zenodo.4613881)   
  - [models.tar.gz](https://zenodo.org/record/4613881/files/models.tar.gz?download=1): Models of protein mutants used for the predictions.  
  - [predictions.tar.gz](https://zenodo.org/record/4613881/files/predictions.tar.gz?download=1): ΔΔG predictions by using DDGun and DDGun3D methods.  
  


#### INSTALLATION

     # Automatic Installation
       git clone https://github.com/biofold/ddgun
       cd ddgun
       python setup.py

     # Manual Installation
     1) Download DDGun
        git clone https://github.com/biofold/ddgun

     2) Install hhblits
        cd ddgun/utils
        git clone https://github.com/soedinglab/hh-suite.git
        mkdir -p hh-suite/build && cd hh-suite/build
        cmake -DCMAKE_INSTALL_PREFIX=.. ..
        make -j 4 && make install

     3) Download uniclust30_2018_08_hhsuite (~25Gb)
        cd ../../../data
        wget http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz

     4) untar uniclust30_2018_08_hhsuite.tar.gz
        tar -xzvf uniclust30_2018_08_hhsuite.tar.gz
        cd ../


#### USAGE

    - Run DDGun 3D:
        ./ddgun_3d.py test/1aar.pdb A test/1aar.muts

        #PDBFILE        CHAIN   VARIANT S_DDG   T_DDG
        1aar.pdb        A       K6A     0.5     0.5
        1aar.pdb        A       T7A     -0.3    -0.3
        1aar.pdb        A       L8A     -0.4    -0.4
        1aar.pdb        A       T9A,G10A        -0.1,-0.1       -0.1


    - Run DDGun Seq:
        ./ddgun_seq.py test/1aar.pdb.A.fasta test/1aar.muts

        #SEQFILE        VARIANT S_DDG   T_DDG
        1aar.pdb.A.fasta        K6A     0.4     0.4
        1aar.pdb.A.fasta        T7A     -0.5    -0.5
        1aar.pdb.A.fasta        L8A     -1.2    -1.2
        1aar.pdb.A.fasta        T9A,G10A        -0.7,-0.6       -0.7


    - Output Legend:
        PDBFILE: Inoput PDB File.
        CHAIN:   Input PDB File Chain.
        SEQFILE: Input Sequence File.
        VARIANT: Comma-separated protein variant in the format XPOSY 
                 (X=Wild-Type Residue, POS=Position, Y=Mutant Residue).
        S_DDG:   Comma-separated predicted DDG of unfolding for single mutants.
        T_DDG:   Final predicted DDG of unfolding. For multiple mutant is 
                 obtained as a combination of the single mutant predictions.

