# Kernels4graphs

## Description

Kernels4graphs implements several kernel based methods for classification of graphs. A basic SVM algorithm and 4 different graph Kernels implementations are provided. The different approaches have been applied to Structure-Activity Relationship Analysis (SAR) to predict the molecular activity of differente chemical compounds described as labeled graphs. The different approaches have been compared and the results are presened in the file report.pdf.
The kernels implemented here are:

+ Nth order walk kernel 
+ Geometric walk kernel 
+ Markovian random walk kernel
+ Subtree kernel


## Dependencies

+ Networkx
+ Numpy
+ sklearn
+ rdkit
+ cvxopt

## References

+ Ralaivola, L., Swamidass, S. J., Saigo, H., & Baldi, P. (2005). Graph kernels for chemical informatics. Neural networks, 18(8), 1093-1110.
+ Helma, C., Cramer, T., Kramer, S., & De Raedt, L. (2004). Data mining and machine learning techniques for the identification of mutagenicity inducing substructures and structure activity relationships of noncongeneric compounds. Journal of chemical information and computer sciences, 44(4), 1402-1411.
+ Kashima, H., Tsuda, K., & Inokuchi, A. (2004). Kernels for graphs. Kernel methods in computational biology, 39(1), 101-113.
+ Mahé, P., Ueda, N., Akutsu, T., Perret, J. L., & Vert, J. P. (2005). Graph kernels for molecular structure-activity relationship analysis with support vector machines. Journal of chemical information and modeling, 45(4), 939-951.
+ Shervashidze, N., & Borgwardt, K. M. (2009). Fast subtree kernels on graphs. In Advances in Neural Information Processing Systems (pp. 1660-1668).
+ Mahé, P., & Vert, J. P. (2009). Graph kernels based on tree patterns for molecules. Machine learning, 75(1), 3-35.