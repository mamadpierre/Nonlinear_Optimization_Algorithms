function [x,f,f_0,g,H,d_0,counters,norms,i,out_null,out_line]= initials(p,x,algo,i)

% function [x,f,f_0,g,H,d_0,counters,norms,i,out_null,out_line]= initials(p,x,algo,i)
%
% Author      : Mohammad Pirhooshyaran
% Description : initial inputs for package
% Input       : p ~ problem
%               i ~ input structure
%               algo ~ one of eight algorithms
%               x ~ initial iterate
% Output      : x ~ final iterate
% counters    : different counting within a structure
% counters    : different norms within a structure


% Termination check
if ~isfield(i,'maxiter')
    i.maxiter = 1e+3; % iteration limit
end

if ~isfield(i,'opttol')
    i.opttol = 1e-06; % optimality toleration
end

%/////////////////////////////
% Line search paramters
if ~isfield(i,'c1ls')
    i.c1ls = 0.0001;
end

if ~isfield(i,'c2ls')
    i.c2ls = 0.9;
end

if ~isfield(i,'amax')
    i.amax = 10;
end

% Backtrack paramter
if ~isfield(i,'rho')
    i.rho = 0.5;
end
%////////////////////////////

%////////////////////////////
% Trust region parameters
if ~isfield(i,'c1tr')
    i.c1tr = 0.1;
end

if ~isfield(i,'c2tr')
    i.c2tr = 0.75;
end

if ~isfield(i,'delta_increase')
    i.delta_increase = 1.5;
end

if ~isfield(i,'delta_decrease')
    i.delta_decrease = 0.5;
end
%//////////////////////////////////

%optimality tolerance for cg
if ~isfield(i,'cgopttol')
    i.cgopttol = 1e-06;
end

%max iteration for cg
if ~isfield(i,'cgmaxiter')
    i.cgmaxiter = 1e+03;
end

% Trust region with SR1 updates
if ~isfield(i,'sr1updatetol')
    i.sr1updatetol = 1e-06;
end

% BFGS Hessian approximation update parameter
if ~isfield(i,'bfgsdamptol')
    i.bfgsdamptol = 0.2;
end

%initializing variables, vectors and matrices needed for starting
%different methods


% Store output strings (Here different outputs are printed for linesearch methods versus trust region methods)
out_line = '==================================================================';

if strcmp(algo,'steepestbacktrack')==1 || strcmp(algo,'steepestwolfe')==1 ...
        || strcmp(algo,'newtonbacktrack')==1 || strcmp(algo,'newtonwolfe')==1 ...
        || strcmp(algo,'bfgsbacktrack')==1 || strcmp(algo,'bfgswolfe')==1
    
    out_data = '  k             f             ||g||             ||d||            D = g^T*d         alpha      #LS_solved        #gradient_evaluation';
    out_null =                                   '----------  -----------  ----------   ----------   ----------';
elseif strcmp(algo,'trustregioncg')==1 || strcmp(algo,'sr1trustregioncg')==1
    out_data = '  k            f             ||g||         ||d||         D = g^T*d        delta       CG iter   #LS_solved    #gradient_evaluation';
    out_null =                                  '----------  -----------  ----------  ----------   ----------   ----------';
end
% Print output header
fprintf('%s\n%s\n%s\n',out_line,out_data,out_line);

% Initialize iteration counter
counters.k = 0;
% Initialize CG subproblem iteration counter
counters.CG_k = 0;
% Initialize counter for number of function evaluation
counters.f = 0;
% Initialize counter for number of gradient evaluation
counters.g = 0;
% Initialize counter for number of hessian evaluation
counters.H = 0;
% Initialize counter for number of linear systems evaluation
counters.LS = 0;
% Initialize SR1 updates counter
counters.SR1 = 0;

% Evaluate function value at x
f = feval(p,x,0);
counters.f =counters.f +1;

% Store initial function value
f_0 = f;

% Evaluate gradient value at x
g = feval(p,x,1);
counters.g =counters.g +1;

% Evaluate gradient norm at x
norms.g = norm(g);

% Store initial gradient norm
norms.g_0 = norms.g;
size_gradient= size(g,1);

% setting two things for trust region methods:
%/ Initialize the search direction for trustregion methods
%// Initialize first trust region radius for trustregion methods

%put the first trust region radius
i.delta = 0.01*norms.g;

% it can also be produced randomly
%(i.deltamax-i.deltamin)*rand(1,1) + i.deltamin;

%Initialize the search direction for trustregion methods
d_0 = zeros(size_gradient,1);

%if strcmp(algo,'trustregioncg')==1 || strcmp(algo,'sr1trustregioncg')==1
%d_0 = zeros(size_gradient,1);
%end

% Initialize hessian matrix for bfgs methods
%if strcmp(algo,'bfgsbacktrack')==1 || strcmp(algo,'bfgswolfe')==1

%start with a coefficient of identity matrix
%H = eye(size_gradient);
%counters.H =counters.H +1;
%end
%this line is just for running the function with output H
H = eye(size_gradient);
%  if  strcmp(algo,'sr1trustregioncg')==1
%  % Evaluate Hessian of x
%         H = feval(p,x,2);
%  end
end


