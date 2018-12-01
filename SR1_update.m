function [x,f,g,H,counters,norms,D,d,i]= SR1_update(p,x_n,x,f,g,H,counters,norms,i,D,d)

% Author      : Mohammad Pirhooshyaran
% Description : specific function for SR1Update
% (reject/accept function)
% Input       : p ~ problem
%               i ~ input structure
%               x ~ initial iterate
% Output      : x ~ final iterate
% counters    : different counting within a structure
% counters    : different norms within a structure


% Evaluate the gradient of x_k+d_k
g_n = feval(p,x_n,1);
counters.g =counters.g +1;
% Find the difference between the step x_k and x_k+1
s = x_n - x;
%s=s';
% Find the difference between the gradients of x_k and x_k+1
y = g_n-g;

% SR1 update

% If we do not accept d, we must skip the SR1update
if (norm(s)~= 0)
    
    %variables needed for sr1 update
    H_S = H*s;
    if (abs((y - H_S)'*s) >= (i.sr1updatetol* norm(y-H_S)*norm(s)))
        
        % Update the Hessian
        H = H + ((y-H_S) * (y-H_S)')/((y-H_S)'*s);
        
        %  else Skip the update in Hessian H = H;
    end
end

x = x_n;

% Update the objective gradient
g = g_n;

% Evaluate gradient norm
norms.g = norm(g);
end