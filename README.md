# Nonlinear_Optimization_Algorithms
A MATLAB Package for Nonlinear Optimization Algorithms

This package contains basic MATLAB implementations of:
1) Steepest Descent algorithm with backtracking
2) Steepest Descent algorithm with Wolf-condition
3) Newton algorithm with backtracking
4) Newton algorithm with Wolf-condition
5) BFGS algorithm with backtracking
6) BFGS algorithm with Wolf-condition
7) Trust Region with Conjugate Gradient (CG) as Subproblem Solver
8) SR1 update Trust Region with Conjugate Gradient (CG) as Subproblem Solver

Implementations are based on: Nocedal, J., & Wright, S. J. (2006). Numerical optimization 2nd.

To run the opt.m optimizer you need to execute the following command in the directory of package:
```
output = opt('problem_name',starting_point,'algorithm',parameters)
```
for instance:
```
[x]= opt('rosenbrock',[1.8;3],'sr1trustregioncg',struct)
[x]= opt('genhumps',[2;2;2;2;2],'newtonbacktrack',struct)
```
If you want to execute all the algorithms for a given problem and starting point, you may use runner.m:
```
runner('problem_name',starting_point)
```

for instance:
```
runner('leastsquares',[3;3;3;3])
```
##additional details:

1) If you want to add new problems of your own, you need to make sure the 'opt' can access function evaluation (case 0), gradient evaluation (case 1) and if needed by the algorithm Hessian evaluation (case 2). 

2) If you wish to change the parameters defualt values, you may add them to structure before running the opt.m
