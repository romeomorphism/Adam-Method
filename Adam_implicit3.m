function [x, y] = Adam_implicit3(f_ode, xRange, yInitial, numSteps)
% Adam third-order implicit method to solve ode dy/dx = f_ode(x, y)
% xRange = [x1, x2] where the solution is approximated in interval [x1, x2]
% yInitial = column vector of initial values for y at x1
% numSteps = number of equal stepsizes from x1 to x2
% x = row vector of value of x whose k-th elments is x1 + (k-1)*stepsize
% y = matrix whose k-th column is the approximate solution at x(k)

% Reminder: y_{n+1} = y_{n} + h*(5/12*f_{n+1} + 8/12*f_n - 1/12f_{n-1})
% For the Cauchy problem dy/dx = -y - 3 * x with y(0) = 1, we need to solve
% a linear equation which is shown in the code
x(1) = xRange(1);
h = (xRange(2) - xRange(1)) / numSteps; %step size
y(:,1) = yInitial;

k = 1;
fValue = f_ode(x(k),y(k));
% Set the exact solution to the initial value
x(k+1) = x(k) + h;
y(:, k+1) = exmp_sol(x(k+1));

for k = 2:numSteps
    fValueold = fValue;
    fValue = f_ode(x(k), y(:,k));
    x(k+1) = x(k) + h;
    y(:,k+1) = y(:,k)/(1+5/12*h) + h/(1+5/12*h)*(-5/4*x(k+1) + 2/3*fValue - 1/12*fValueold);
end

