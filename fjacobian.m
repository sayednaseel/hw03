function dfdx = fjacobian(f, x0, dx)
% FJACOBIAN   Computes Jacobian (partial derivative) of function
%
% dfdx = fjacobian(@f, x0 [, dx])
%
%   Uses finite differences to compute partial derivative of vector 
%   function f evaluated at vector x0. 
%   Optional argument dx specifies finite difference
%   step. Note that argument f should typically be entered as a
%   function handle, e.g. @f.
%
%   When dx is not given or empty, a default of 1e-6 is used.

  if nargin < 3 || isempty(dx)
    dx = 1e-6;
  end
  
  f0 = f(x0);
  J = zeros(length(f0), length(x0));
  for i = 1:length(x0)
    xperturbed = x0;
    xperturbed(i) = xperturbed(i) + dx;
    df(:,i) = f(xperturbed) - f0;
  end
  dfdx = df / dx;
end