function     [x,f,g,H,counters,norms,i]=  stepacceptance(p,x,algo,f,g,H,counters,norms,i,D,d)


% Author      : Mohammad Pirhooshyaran
% Description : step size fucntion for different algorithms
% Input       : p ~ problem
%               i ~ input structure
%               algo ~ one of eight algorithms
%               x ~ initial iterate
% Output      : x ~ final iterate
% counters    : different counting within a structure
% counters    : different norms within a structure





% Run line search for linesearch methods and trust region for trust region
% methods

% put the first step length
a=1;

%different line search algorithms
if strcmp(algo,'steepestbacktrack')==1 || strcmp(algo,'steepestwolfe')==1 ...
        || strcmp(algo,'newtonbacktrack')==1 || strcmp(algo,'newtonwolfe')==1 ...
        || strcmp(algo,'bfgsbacktrack')==1 || strcmp(algo,'bfgswolfe')==1
    
    %line search
    [a,f,counters,i] = linesearch(p,x,f,d,D,algo,counters,i);
    %print the steplenght
    fprintf('%.4e\t\t\t',a);
    
    
    %specific update for just two bfgs
    if strcmp(algo,'bfgsbacktrack')==1 || strcmp(algo,'bfgswolfe')==1
        
        %updating function for two bfgs
        [x,f,a,g,H,counters,norms,D,d,i]= bfgs_update(p,x,f,g,H,counters,norms,i,D,d,a);
        
        %other linesearch algorithms
    elseif strcmp(algo,'steepestbacktrack')==1 || strcmp(algo,'steepestwolfe')==1 ...
            || strcmp(algo,'newtonbacktrack')==1 || strcmp(algo,'newtonwolfe')==1 ...
            
        %update for other linesearch algorithms
        
        % Update iterate
        x = x + a*d;
        
        % Evaluate objective gradient
        g = feval(p,x,1);
        counters.g = counters.g + 1;
        
        % Evaluate gradient norm
        norms.g = norm(g);
        
    end
    
    
    %different trustregion algorithms
elseif strcmp(algo,'trustregioncg')==1 || strcmp(algo,'sr1trustregioncg')==1
    
    %trustregion function
    [x_n,x,f,d,D,g,H,counters,i] = trustregion(p,x,f,d,D,g,H,counters,i);
    
    if strcmp(algo,'trustregioncg')==1
        
        x=x_n;
        
        % Evaluate objective gradient
        g = feval(p,x,1);
        counters.g = counters.g + 1;
        
        % Evaluate gradient norm
        norms.g = norm(g);
        
    elseif strcmp(algo,'sr1trustregioncg')==1
        
        %Specific update for sr1
        [x,f,g,H,counters,norms,D,d,i]= SR1_update(p,x_n,x,f,g,H,counters,norms,i,D,d);
    end
    
    
    
    % Print trust region radius and CG iteration
    fprintf('%.4e\t\t%d\t\t',i.delta, counters.CG_k);
end

% Print trust region radius and CG iteration
fprintf('%d\t\t\t %d \n',counters.LS, counters.g);
end
