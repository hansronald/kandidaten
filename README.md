# kandidaten

Code for bachelor thesis High-dimensional inverse covariance estimation - evaluation of new method using Graphical Lasso.

Abstract
Partial correlation between variables can be obtained implicitly through the precision matrix.
To estimate the precision matrix, the covariance matrix, computed from data, can be inverted.
Problems arise when the number of variables p is larger than the number of observations $n$,
since the covariance matrix then becomes low-rank and can not be inverted.
Former methods to solve this problem are uncertain and therefore a new method called \textit{k-glasso} has been developed.
In October 2017, the article  "Improving the Graphical Lasso Estimation for the Precision Matrix Through Roots of the Sample Covariance Matrix" (Avagyan 2017)
was published, where the method k-glasso was presented. The claim was that the method outperforms previous methods.
The aim of this study was to examine the performance of k-glasso on large,
sparse precision matrices that are similar to networks from applications. For the simulation,
two block diagonal precision matrices were generated for different values of p, where the underlying network had a
scale-free distribution. Two variations of one of the models from the original article was also investigated.
The k:th root of the empirical covariance matrix was calculated by taking the k:th root of the diagonal in its
eigenvalue decomposition. The R-function huge() was used to calculate the estimates of k-glasso for a certain
sequence of lambda. Then the data was transformed back by taking the estimate to the power of k.
By doing 100 replicates, the mean of the evaluation measures were computed.
The method was also applied to cancer data from real observations. The results in this study were not consistent
with the results from the original article since the performance of k-glasso seems to be data dependent.