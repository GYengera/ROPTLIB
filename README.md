# ROPTLIB with SE-Sync Optimization Example

The main function for the SE-Sync optimization problem is in /test/TestStieSync.cpp
The computation of the function, gradient and hessian can be found in /Problems/StieSync/ folder.

The above mentioned code is self-explanatory as all necessary comments are included.

To compile this example use the command: $ make ROPTLIB TP=TestStieSync

After compilation, the command to run the example is the following:
$ ./TestStieSync <Q Dataset Path> <n> <d> <p> <Initial value X0 path>

If the path is not provided to the data for initializing the value of X, then a random point is created on the stacked stiefel manifold.

The optimum value of X is stored in a text file. The path to the file can be altered within the code.
Gradient and hessian checks are displayed in the terminal for verification.
