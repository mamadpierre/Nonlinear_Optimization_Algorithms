function [a,f,counters,i] = linesearch(p,x,f,d,D,algo,counters,i)

% function [a,f] = linesearch(p,x,f,d,D,o)
%
% Revised by  : Mohammad Pirhooshyaran
% Description : Line search subroutines for steepestdescent.m.
% Input       : p ~ problem function handle
%               x ~ point
%               f ~ function value
%               d ~ search direction
%               D ~ directional derivative value
%               o ~ line search option
%               i ~ input parameters
%               counters ~ counter of the function, gradient and hessian
% Output      : a ~ step size
%               f ~ updated function value

if strcmp(algo,'steepestbacktrack')==1 || strcmp(algo,'newtonbacktrack')==1 || strcmp(algo,'bfgsbacktrack')==1
    
    % put the first step length
    a = 1;
    
    % Initialize iteration counter
    iter = 0;
    
    %Iteration main loop
    while 1
        
        % put objective value for iteration
        f_n = feval(p,x+a*d,0);
        counters.f =counters.f +1;
        
        % Check Sufficient condition (Armijo condtion)
        if(f_n > (f + a*D*i.c1ls))
            
            % Reduction of steplength till satisfaction of sufficient decrease
            a = i.rho*a;
            
        else
            % put the objective value when it satisfies the Armijo condition
            f = f_n;
            break
            
        end
        
        % Increment iteration counter
        iter = iter+1;
        
    end
    
    % Wolf line search
elseif strcmp(algo,'steepestwolfe')==1 || strcmp(algo,'newtonwolfe')==1 || strcmp(algo,'bfgswolfe')==1
    
    % Set max and min step lengths
    i.amax = 10;
    a_0 = 0;
    
    % Set the first steplength
    a = 1; %(i.amax-a0)*rand(1,1) + a0;
    
    
    % Initialize iteration counter
    iter = 1;
    
    %Iteration main loop
    while 1
        
        % Evaluate the function value of the iterate
        f_n = feval(p,x+a*d,0);
        counters.f =counters.f +1;
        
        % Check sufficient decrease condition
        if (f_n > f + i.c1ls*a*D) || ((iter > 1) && (f_n >= feval(p,x+a_0*d,0)))
            counters.f =counters.f +1;
            
            % Run zoom function and terminate the code with the step length
            [a_opt, counters] = zoom(a_0,a,f,i.c1ls,i.c2ls,D,p,x,d,counters);
            break
            
        end
        
        % Evaluate the gradient of the iterate
        g_n = feval(p,x+a*d,1);
        counters.g =counters.g +1;
        
        % Check curvature condition (here we used strong curviture
        % condition)
        if (abs(g_n'*d) <= -i.c2ls * D)
            
            % Find the optimum steplength then stop
            a_opt = a;
            break
            
        end
        
        % Check the last condition
        if (g_n'*d >= 0)
            
            % Run zoom function and terminate the code with the step length
            [a_opt, counters] = zoom(a,a_0,f,i.c1ls,i.c2ls,D,p,x,d,counters);
            break
            
        end
        
        % Update the previous iterate for the next loop
        a_0 = a;
        
        % Choose the new steplength between the maximum and previous
        % iterate
        a = (i.amax+a_0)/2;
        
        % a = (am-a0)*rand(1,1) + a0; % it can also be chosen randomized
        
        %  counter
        iter = iter+1;
        
    end
    
    % find the final steplength and function value
    a = a_opt;
    f = feval(p,x+a*d,0);
    counters.f =counters.f +1;
    
end
end

% ZOOM FUNCTION checks three prorerties below:
%1) sufficient decrease property
%2) strong curviture property
%3) derivitive of steplengh is greater than or equal to zero
%=================================================
function [a_star, counters] = zoom(a_L,a_U,f,c1ls,c2ls,D,p,x,d,counters)

% function a = zoom()
% Description : Zoom function for Wolfe line search.
% Input       : a_L ~ lower bound for the step length
%               a_U ~ higher bound for the step length
%               x ~ point
%               f ~ function value
%               d ~ search direction
%               D ~ directional derivative value
%               p ~ problem function handle
% Output      : a_star ~ step length



% Initialize iteration counter
iter_zoom = 0;

%Iteration loop
while 1
    
    % Calculate the trial step length between the bounds
    a = (a_U+a_L)/2;
    % a = (ah-al)*rand(1,1) + al; % it can also be chosen randomized
    
    % Evaluate the function value of the iterate
    f_n =  feval(p,x+a*d,0);
    counters.f =counters.f +1;
    
    % Check sufficient decrease condition (Property (1))
    if  (f_n > f + a*D*c1ls) || (f_n >= feval(p,x+a_L*d,0))
        
        counters.f =counters.f +1;
        
        % find the higher bound
        a_U = a;
        
    else
        
        % Evaluate the gradient of the iterate
        g = feval(p,x+a*d,1);
        counters.g =counters.g +1;
        
        % Check the curvature condition (Property (2))
        if (abs(g'*d) <= -c2ls * D)
            
            % Find the steplength satisfying the properties and stop
            a_star = a;
            break
            
        end
        
        % Check Property (3)
        if (g'*d*(a_U - a_L) >= 0)
            
            % Update the higher bound
            a_U = a_L;
            
        end
        
        % Update the lower bound
        a_L = a;
        
    end
    iter_zoom = iter_zoom+1;
    
end
end








