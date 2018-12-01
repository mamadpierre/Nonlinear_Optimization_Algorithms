function [x,f,a,g,H,counters,norms,D,d,i]= bfgs_update(p,x,f,g,H,counters,norms,i,D,d,a)

% Author      : Mohammad Pirhooshyaran
% Description : specific function for BFGS update
% (reject/accept function)
% Input       : p ~ problem 
%               i ~ input structure
%               x ~ initial iterate
% Output      : x ~ final iterate
% counters    : different counting within a structure
% counters    : different norms within a structure



% Evaluate the gradient of the next itrate
    g_n = feval(p,x+a*d,1);
    counters.g =counters.g +1;
    
    % evaluate the difference between the step x_k and x_k+1
    s = a*d;
    
    % evaluate difference between the gradients of x_k and x_k+1
    y = g_n-g;
    
    % Define the scalar theta
    if ((s'*y) >= i.bfgsdamptol *s'*H*s)
        
      theta = 1;
      
    else
        
      theta = ((1-i.bfgsdamptol)*s'*H*s)/(s'*H*s - s'*y);
      
    end

    % substitue y_k with r_k
    r = theta * y + (1-theta)*H*s;

    % Update and substitue the Hessian of the next iterate
    H = H -( H*s*s'*H)/(s'*H*s)+(r*r')/(s'*r);

    % Update and substitute the gradient of the next iterate amd its norm
    g = g_n;
    norms.g = norm(g);
    
    % find the next iterate
    x =  x + a*d;
end