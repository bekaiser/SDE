function [ k ] = sde_RK4(T,z_dt,A,G)
% RK4 k coefficient generator. Uses the algorithm
% described by Hansen & Penland (2006).

% Here, the spatial correlation matrix G is 
% integrated for uncorrelated, autonomous
% spatial correlations. The dynamical operator
% A is also autonomous, so A*T is autonomous.

% input variables:
% size(T) = [1 N], N = number of state variables
% t = a scalar variable; the time
% dt = the time step.

k = A*T+G*z_dt;

end
