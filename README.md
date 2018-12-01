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

To run the 'opt' optimizer you need to execute the following command in the directory of package:
```
output = opt('problem_name',starting_point,'algorithm',parameters)
```
for instance:
```
[x]= opt('rosenbrock',[1.8;3],'sr1trustregioncg',struct)
[x]= opt('genhumps',[2;2;2;2;2],'newtonbacktrack',struct)
```
If you want to execute all the algorithms for a given problem and starting point, you may use runner.m:

for instance:
```
runner('leastsquares',[3;3;3;3])
```


