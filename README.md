# Out-of-Sample-Extension

# example1.m

Dowload the poker dataset and put it in the same directory of this .m file.
To see the plots of chapter 3, you may use the following command:

example1('poker', 20000, [50:50:800])

A total of three graphs should be plotted:
1. size of difference matrix vs. k
2. 1 - largest eigenvalue of difference matrix vs. k
3. 5 other eigenvalues vs. k

Also, a table consisitng of the eigevalues of the difference matrix for different values of k will be printed

If you wish to use your own dataset, you may use this as follows:

example1(your_own_dataset_as_matlab_matrix, number_of_samples_to_use, vector_consisitng_of_k_values)

note that the data matrix is in the format samples * features.


#############

# example2.m
