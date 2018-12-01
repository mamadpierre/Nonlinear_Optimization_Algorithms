function [x_n,k] = cg(delta,A,b,x)

% function xn = cg(delta,A,b,x)
%
% Description : Conjugate Gradient method.
% Input       : delta ~ trust region radius
%               A ~ symmetric matrix
%               b ~ right-hand side vector
%               x ~ initial iterate
% Output      : xn ~ final iterate
% Revised by  : MP

% Evaluate initial residual
r = A*x-b;

% Evaluate residual norm
norms.r0 = norm(r);

% Evaluate initial direction
p = -r;

% Evaluate direction norm
norms.p = norm(p);

% Initialize iteration counter
k = 0;

% Main LOOP of CG
while 1
    
    % calculate matrix-vector product
    Ap = A*p;
    
    % calculate vector-vector product
    pAp = p'*Ap;
    
    if (pAp < 0)
        pp=p'*p;
        % Here we evaluate the positive value of alpha in a way that ||x_k+alpha_k*p_k||=delta_k
        alpha = max(roots([(pp) (2*x'*p) (x'*x- delta.^2)]));
        
        % find the next iterate
        x_n = x+ alpha*p;
        break
        
    else
        
        alpha = (r'*r)/pAp;
        
    end
    
    if norm(x+alpha*p) > delta
                pp=p'*p;
        % Here we evaluate the positive value of alpha in a way that ||x_k+alpha_k*p_k||=delta_k
        alpha = max(roots([(pp) (2*x'*p) (x'*x- delta.^2)]));
        
        % find the next iterate
        x_n = x+ alpha*p;
        break
        
    else
        
        % find the next iterate and residual
        x_n = x + alpha*p;
        r_n = r + alpha'*Ap;
        
    end
    
    % when next residual is close enough to zero i.e. ||r_k+1||~0
    if norm(r_n)<= 0.000001*max(norms.r0,1)
        
        % update the next iterate and stop
        x_n = x + alpha*p;
        break
        
    else
        
        % find the initial direction and beta
        beta = (r_n'*r_n)/(r'*r);
        p = -r_n+ beta*p;
        
    end
    
    % find the next iterate and residual for updating new loop
    r=r_n;
    x=x_n;
    
    %  counter plus
    k = k + 1;
    
end

end


