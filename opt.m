function [x]= opt(p,x,algo,i)

% function package [x]= opt(p,x,algo,i)
%
% Author      : Mohammad Pirhooshyaran
% Description : package of eight different algorithms
% Input       : p ~ problem like 'rosenbrock'
%               i ~ input structure
%               algo ~ one of eight algorithms 'sr1trustregioncg'
%               x ~ initial iterate
% Output      : x ~ final iterate



tic

%initialize
[x,f,f_0,g,H,d_0,counters,norms,i,out_null,out_line]= initials(p,x,algo,i);

%main loop iteration
while 1
    % Check termination conditions
    if counters.k > i.maxiter || norms.g <= i.opttol*max(norms.g_0,1), break; end;
    
    %evaluate decrease direction and directional derivitives
    [x,f,g,H,counters,norms,D,d,i]=  searchdirection(p,x,algo,f,g,H,counters,norms,i,d_0);
    
    % for different trust region algorithms
    if strcmp(algo,'trustregioncg')==1 || strcmp(algo,'sr1trustregioncg')==1
        
        %check CG iteration limit and optimality tolerance for two
        if counters.CG_k > i.cgmaxiter || norms.d <= i.cgopttol, break; end;
    end
    
    % acceptance and update function (here a linesearch or trust region
    % will implement)
    [x,f,g,H,counters,norms,i]=  stepacceptance(p,x,algo,f,g,H,counters,norms,i,D,d);
    
    %update iterate
    counters.k = counters.k + 1;
    
end

% Print output footer
fprintf('%s\n%s\n',out_null,out_line);
result= zeros(1,4);
result(1,1)=counters.k-1;
result(1,2)=f;
result(1,3)=norms.g;
result(1,4)=toc;
%result
%more results
if (norms.g > i.opttol) && (counters.k > i.maxiter)
    fprintf('%s%s%s%s%s%d\n','During solving problem ', p,' with algorithm ' ,algo,' We reached the maximum iteration allowed which is ', i.maxiter)
else
    fprintf('%s%s%s%s%s%d%s\n','The problem ', p,' was solved with algorithm ' ,algo,' within ', counters.k-1,' iterations')
end

fprintf('%s\n','The algorithm reached possible stationary point')
fprintf('%12.8f\n' ,x);
fprintf('%s%12.8f%s%12.8f\n','Algorithm starts with initial objective value ' ,f_0 ,' and finishes with final objective value ', f);
fprintf('%s%12.8f%s%12.8f\n','Algorithm starts with initial gradient norm value ' ,norms.g_0 ,' and finishes with final gradient norm value ', norms.g);
end



