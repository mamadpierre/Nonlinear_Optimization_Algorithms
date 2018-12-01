function runner(problem,x)
   


% Set input parameters

    % Termination check
    i.maxiter = 1e+3; % iteration limit
    i.opttol = 1e-06; % optimality condition
    
    % Line search paramters
    i.c1ls = 0.0001; %Armijo parameter
    i.c2ls = 0.9; % Strong curvature condition parameter
    i.amax = 10; % max step size
    
    % Backtrack paramter
    i.rho = 0.5; 
    
    % TR parameters
    i.c1tr = 0.1; 
    i.c2tr = 0.75;
    i.delta_increase = 1.5; % increae amount
    i.delta_decrease = 0.5; % decrease amount

    
    % Cojugate Gradient parameters
    i.cgopttol = 1e-06; % Optimality tolerance
    i.cgmaxiter = 1e+03; % iteration limit for CG
    
    %  SR1 updates
    i.sr1updatetol = 1e-06; % parameter SR1 hessian approximation
    
    % BFGS HEssian approximation updates
    i.bfgsdamptol = 0.2; % tolerance for BFGS hessian approximation
    
    % Run steepest descent with backtracking line search
    x_n = opt(problem, x, 'steepestbacktrack',i);
    
    % Run steepest descent with Wolfe line search
    x_n = opt(problem, x, 'steepestwolfe',i);
    
    % Run Newton with backtracking line search
    x_n = opt(problem, x, 'newtonbacktrack',i);  
    
    % Run Newton with Wolfe line search
    x_n = opt(problem, x, 'newtonwolfe',i);   
    
    % Run BFGS with backtracking line search
    x_n = opt(problem, x, 'bfgsbacktrack',i);
    
    % Run BFGS with Wolfe line search
    x_n = opt(problem, x, 'bfgswolfe',i);
    
    % Run trust region method with conjugate gradient subproblem
    x_n = opt(problem, x, 'trustregioncg',i);
    
    % Run SR1 updated trust region method with conjugate gradient subproblem
    x_n = opt(problem, x, 'sr1trustregioncg',i);
    

    
end