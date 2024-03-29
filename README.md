### Computing 2pCF and MCF

CODE : wp_mp_compute.c

This code computes projected two-point and marked correlation functions of real data and jackknife (JK) samples. It reads real galaxies from the file _real\_galaxies_ and random galaxies from the file _random\_galaxies_. 

Using random galaxies, it computes RR and saves as _RR.txt_ in the automatically generated directory _biproducts_. If a _RR.txt_ file exists already, that will be used always. The same happens for jackknife RR files as well.

Using real galaxies along with their properties, the code computes the MCF along with the 2pCF. The MCF computation for each property is done parallely. The weights of random galaxies given in the input file is made use of. 

Once real data is done, computations of JK are done parallely. The real and random JK data files are read from _jackknife\_data_ directory. The JK files should be of name _jk\<jknr\>\_real\_galaxies_ and _jk\<jknr\>\_random\_galaxies_.

A _result_ directory will be generated with many files. The final result file of interest would be _wpRealAll.txt_, _mpRealAllShuffle\_mk\<nr\>.txt_ and _mpRealAll\_mk\<nr\>.txt_. _wpRealAll.txt_ and _mpRealAll\_mk\<nr\>.txt_ contains CFs in real data and its JK samples.  _mpRealAllShuffle\_mk\<nr\>.txt_ contains MCF of real sample with the 100 shuffles. 

The formats of all files are given below.

#### REQUIRED FILES/FOLDERS:

1. Header file: _my\_functions.h_ 
2. Real data file : _real\_galaxies_ 
3. Random data file : _random\_galaxies_
4. Directory _jackknife\_data_ with real and random JK files
5. List of properties: _propertylist.txt_ . This file should have the headers of property columns in the real data file.

#### FILE COLUMNS (seperated by tabspace, all input files should have the header):

1. All real and JK real data files: #galaxy\_id, ra, dec, redshift, prop1, prop2, ...
2. All random and JK random data files : galaxy\_id, ra, dec, redshift, weight 
3. Result file wpRealAll.txt : rp, wp\_real, wp\_copy1, wp\_copy2, ... #copy: jk
  Result file mpRealAll_mk\<nr\>.txt : rp, Mp_real, Mp_copy1, Mp_copy2, ... #copy: jk
4. Results file mpRealAllShuffle_mk\<nr\>.txt : rp, Mp_real, Mp_shuffle1, Mp_shuffle2, ...

#### PARAMETERS TO BE SET:
 
- non_prop : Number of columns in real file that are not to be used as marks. These columns should be the first columns.
- n_jacks : number of JK samples present in the folder 'jackknife_data' (generates error if both n_boots and n_jacks have non-zero values)
- n_shuffle : number of times the marks are shuffled to find error in MCF
- bin_rp : number of bins in rp
- bin_pi : number of bins in pi
- rp_init : centre of first (smallest) bin in rp
- pi_init : centre of first (smallest) bin in pi
- dlrp : binwidth in rp (log scale)
- dpi : binwidth in pi (linear scale)

#### TO RUN:

`gcc -fopenmp -o wp_mp_compute wp_mp_compute.c -lm && ./wp_mp_compute`
  
### Fitting 2pCF
  
"fit_code.c" fits the 2pCF with a power-law. But before the fitting procedure, we need to remove NaNs from the result files and do SVD cleaning of the covariance matrix (Sureshkumar+2021). This is done with the python code "svd_invertor.py". In that file, declare if you need to do the SVD cleaning by assigning 1 or 0 for the variable to_filter. Then run the python code `python svd_intertor.py`.
  
After the cleaning, run `gcc -o fit_code fit_code.c -lm && ./fit_code` to fit the CFs and general final results. This action generates a directory _final_ which contails _wp_params.txt_, _final_wp.txt_, _final_Mp_\<propname\>.txt_. 
  
  - _wp_params.txt_ contains the 2pCF power-law parameters r0,error_r0,gamma, error_gamma.
  - _final_wp.txt_ has columns rp, wp(rp), error_wp(rp). 
  - _final_Mp_propname.txt_ has columns rp, Mp(rp), errshuffle_Mp(rp), errJK_Mp(rp)  
  
For information about correlation functions, please read the Chapter 3 of my PhD thesis: https://doi.org/10.5281/zenodo.7572218.

### Sample data to test

In this folder, I provide a sample data - 'real_galaxies' and 'random_galaxies'. Running the code on this sample data should give the results as in Fig. 3.11 (2pCF) and Fig. 3.18 (MCF) of my PhD thesis. See Sect. 3.2 of the thesis to know more about how this dataset was generated from the GAMA survey.

For various applications of 2pCF and MCF, please see the papers:

Sureshkumar et al. 2021: Galaxy and Mass Assembly (GAMA): Tracing galaxy environment using the marked correlation function, A&A 653, A35 (arXiv:2102.04177).

Sureshkumar et al. 2023: Galaxy and Mass Assembly (GAMA): Mid-infrared properties as tracers of galaxy environment, A&A 669, A27 (arXiv:2201.10480).



**For any queries, CONTACT:**

Unnikrishnan Sureshkumar
ukp1513@gmail.com
