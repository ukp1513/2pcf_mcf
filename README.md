CODE : wp_mp_compute.c

This code computes projected two-point and marked correlation functions of real data and bootstrap (BS)/jackknife (JK) samples. It reads real galaxies from the file 'real_galaxies' and random galaxies from the file 'random_galaxies'. 

Using random galaxies, it computes RR and saves as 'RR.txt' in the automatically generated directory 'biproducts'. If a RR.txt file exists already, that will be used always. The same RR file will be used for all BS samples in case of individual galaxy bootstrapping.  

Using real galaxies along with their properties, the code computes the MCF along with the 2pCF. The MCF computation for each property is done parallely. The weights of random galaxies given in the input file is made use of. 

Once real data is done, computations of BS/JK are done parallely. In case of blockwise BS and JK method, the real and random BS/JK data files are read from 'bootstrap_data' and 'jackknife_data' directories respectively. The BS files should be of name 'bs<bsnr>_real_galaxies' and 'bs<bsnr>_random_galaxies'. The JK files should be of name 'jk<jknr>_real_galaxies' and 'jk<jknr>_random_galaxies'.

A 'result' directory will be generated with many files. The final result file of interest would be 'wpRealAll.txt', 'mpRealAllShuffle_mk<nr>.txt' and 'mpRealAll_mk<nr>.txt'. 'wpRealAll.txt' and 'mpRealAll_mk<nr>.txt' contains CFs in real data and its JK/BS samples.  'mpRealAllShuffle_mk<nr>.txt' contains MCF of real sample with the 100 shuffles. 

The formats of all files are given below.

REQUIRED FILES/FOLDERS:

1. Header file: my_functions.h 
2. Real data file : real_galaxies 
3. Random data file : random_galaxies
4. If blockwise bootstrap method, directory 'bootstrap_data' with real and random BS files
5. If jackknife method, directory 'jackknife_data' with real and random JK files

FILE COLUMNS (seperated by tabspace, all input files should have the header):

1. All real and JK/BS real data files: #galaxy_id, ra, dec, redshift, prop1, prop2, ...
2. All random and JK/BS random data files : galaxy_id, ra, dec, redshift, weight 
3. Result file wpRealAll.txt : rp, wp_real, wp_copy1, wp_copy2, ... #copy: bs/jk
  Result file mpRealAll_mk<nr>.txt : rp, Mp_real, Mp_copy1, Mp_copy2, ... #copy: bs/jk
4. Results file mpRealAllShuffle_mk<nr>.txt : rp, Mp_real, Mp_shuffle1, Mp_shuffle2, ...

PARAMETERS TO BE SET:
 
- non_prop : Number of columns in real file that are not to be used as marks. These columns should be the first columns.
- boot_meth : type of bootstrapping (individual or blockwise)
- n_boots : number of BS samples present in the folder 'bootstrap_data'
- n_jacks : number of JK samples present in the folder 'jackknife_data' (generates error if both n_boots and n_jacks have non-zero values)
- n_shuffle : number of times the marks are shuffled to find error in MCF
- bin_rp : number of bins in rp
- bin_pi : number of bins in pi
- rp_init : centre of first (smallest) bin in rp
- pi_init : centre of first (smallest) bin in pi
- dlrp : binwidth in rp (log scale)
- dpi : binwidth in pi (linear scale)

TO RUN:

gcc -fopenmp -o wp_mp_compute wp_mp_compute.c -lm
./wp_mp_compute

For any queries, CONTACT:

Unnikrishnan Sureshkumar
ukp1513@gmail.com
