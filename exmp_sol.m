function [y] = exmp_sol(x)
% This is the exact solution of the ode system we consider in the test
% dy/dx = -y - 3 * x with y(0) = 1
y = -2 * exp(-x) -3 * x + 3;
end

