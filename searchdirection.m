function [x,f,g,H,counters,norms,D,d,i]=  searchdirection(p,x,algo,f,g,H,counters,norms,i,d_0)

%
% Author      : Mohammad Pirhooshyaran
% Description : search direction fucntion for different algorithms
% Input       : p ~ problem
%               i ~ input structure
%               algo ~ one of eight algorithms
%               x ~ initial iterate
% Output      : x ~ final iterate
% counters    : different counting within a structure
% counters    : different norms within a structure
% d           : direction
% D           : directional derivative

% Evaluate search direction (d)
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\%
if strcmp(algo,'steepestbacktrack')==1 || strcmp(algo,'steepestwolfe')==1
    
    % d for steepest descent is just the negative of gradient
    d = -g;
end
%other methods except two bfgs methods needs evaluation of hessian in
%order to find d
%(for bfgs methods we use multiple of identity as H)
if strcmp(algo,'newtonbacktrack')==1 || strcmp(algo,'newtonwolfe')==1 ...
        || strcmp(algo,'trustregioncg')==1 || strcmp(algo,'sr1trustregioncg')==1
    
    % Evaluate Hessian of x
    H = feval(p,x,2);
    counters.H = counters.H + 1;
end

%checking positive definiteness of hessian for 4 methods
if strcmp(algo,'newtonbacktrack')==1 || strcmp(algo,'newtonwolfe')==1 ...
        || strcmp(algo,'bfgsbacktrack')==1 || strcmp(algo,'bfgswolfe')==1
    % calculating the number of variables
    n = length(x);
    
    % initial increment amount
    x_i = 0.0001;
    
    % Approximate Hessian until it becomes positive definite
    while min(eig(H))< 1e-08
        
        % Increment eigenvalues to have positive definite Hessian
        H = H + x_i * eye(n);
        
        % Update the increment by order of ten
        x_i = 10*x_i;
        
    end
end
%direction evaluation for 4 methods
if strcmp(algo,'newtonbacktrack')==1 || strcmp(algo,'newtonwolfe')==1 ...
        || strcmp(algo,'bfgsbacktrack')==1 || strcmp(algo,'bfgswolfe')==1
    
    % Solve linear system
    d = -H\g;
    counters.LS = counters.LS + 1;
    
    %cg running for two trust region methods
elseif strcmp(algo,'trustregioncg')==1 ...
        || strcmp(algo,'sr1trustregioncg')==1
    
    % Run conjugate gradient subproblem
    [d,counters.CG_k] = cg(i.delta,H,-g,d_0);
    
end

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\%


% Evaluate norm of direction
norms.d = norm(d);




% Evaluate directional derivitive  (D)
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\%
% Evaluate directional derivative
D = g'*d;

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\%



% Print iterate information  %Print search direction information
fprintf('%4d\t\t%+.4e\t\t%.4e\t\t',counters.k,f,norms.g);
%fprintf('%4d  %+.4e  %.4e  ',counters.k,f,norms.g);
fprintf('%.4e\t\t%.4e\t\t',norms.d,D);

end

