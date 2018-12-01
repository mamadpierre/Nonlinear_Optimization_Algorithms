function v = rosenbrock(x,o)

% function v = rosenbrock(x,o)
%
% Description : Rosenbrock function evaluator.
% Input       : x ~ current iterate
%               o ~ evaluation option
%                     0 ~ function value
%                     1 ~ gradient value
%                     2 ~ Hessian value
% Output      : v ~ function, gradient, or Hessian value

% Switch on o
switch o
  
  case 0

    % Evaluate function
    v = 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
  
  case 1
  
    % Evaluate gradient
    v = [-x(1)*400*(x(2)-x(1)^2)-2*(1-x(1));
               200*(x(2)-x(1)^2)            ];

  case 2

    % Evaluate Hessian
    v = [-400*(x(2)-3*x(1)^2)+2 -400*x(1);
         -400*x(1)               200      ];
  
end