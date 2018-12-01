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
To run you need to:
'''
output = opt('problem_name',starting_point,'algorithm',parameters)
'''

For convinience there is a runner.m file containing 


