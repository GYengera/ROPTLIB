# SE-Sync Optimization Example with ROPTLIB

SE-Sync is a certifiably correct algorithm for synchronization over the special Euclidean group, a common problem arising in the context of 2D and 3D geometric estimation such as pose-graph SLAM for example. A detailed explanation of the algorithm can be found in this [paper](https://arxiv.org/pdf/1611.00128.pdf).

This code implements the riemannian optimization problem given in problem 4 (equation 13) of the [paper](https://arxiv.org/pdf/1611.00128.pdf). 

The main function for the SE-Sync optimization problem is in /test/TestStieSync.cpp

The computation of the function, gradient and hessian (shown in equation 14 of the [paper](https://arxiv.org/pdf/1611.00128.pdf)) can be found in /Problems/StieSync/ folder.

To compile this example run the following command:

$ make ROPTLIB TP=TestStieSync

After compilation, the command to run the example is the following:

$ < path >/TestStieSync < Q Dataset Path> < n > < d > < p > < Initial value X0 path >

Datasets for matrix Q can be found at: http://mapir.isa.uma.es/jbriales/GSOC_project8_data.zip 

If the path is not provided to the data for initializing the value of X, then a random point is created on the stacked Stiefel manifold.

The optimum value of X is stored in a text file. The path to the file can be altered within the code.
Gradient and hessian checks are displayed in the terminal for verification.
