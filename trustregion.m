function [x_n,x,f,d,D,g,H,counters,i] = trustregion(p,x,f,d,D,g,H,counters,i)


% Author      : Mohammad Pirhooshyaran
% Description : specific function for trust region algorithms
% (reject/accept function)
% Input       : p ~ problem
%               i ~ input structure
%               x ~ initial iterate
% Output      : x ~ final iterate
% counters    : different counting within a structure
% counters    : different norms within a structure


% Evaluation of function value at new point x_k+d_k
f_n = feval(p,x+d,0);
counters.f =counters.f +1;

% Find the predicted reduction in f calculated between x_k and the new point x_k+d_k
m_reduction = -g'*d - 1/2* d'*H*d;

% Find rho which is the actual reduction of f divided by predicted reduction of f
rho = (f - f_n)/m_reduction;

% If  model approximates f accurately and the reduction in f is big
% enough
if (rho >= i.c2tr)
    
    % First accept d_k and update x and function value
    x_n = x + d;
    f = f_n;
    
    % Second inrease the radius
    i.delta = i.delta * i.delta_increase;
    
    % If model may not approximate f accurately but we experience significant reduce in f
elseif (rho < i.c2tr) && (rho > i.c1tr)
    
    % First accept d_k and update x and function value
    x_n = x + d;
    f = f_n;
    
    % Second, maintain the radius
    i.delta = i.delta;
    
    % If we do not have significant reduced in f and model does not
    % approximate f accurately
elseif (rho <= i.c1tr)
    
    % First, reject d_k and keep the function value
    x_n = x;
    f = f;
    
    % Second, reduce the radius
    i.delta = i.delta * i.delta_decrease;
    
end
end