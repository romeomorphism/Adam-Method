function [x, y] = Adam_implicit4(f_ode, xRange, yInitial, numSteps)
% Adam fourth-order implicit method to solve ode dy/dx = f_ode(x, y)
% xRange = [x1, x2] where the solution is approximated in interval [x1, x2]
% yInitial = column vector of initial values for y at x1
% numSteps = number of equal stepsizes from x1 to x2
% x = row vector of value of x whose k-th elments is x1 + (k-1)*stepsize
% y = matrix whose k-th column is the approximate solution at x(k)

% Reminder: y_{n+1} = y_{n} + h*(9/24*f_{n+1} + 19/24*f_n - 5/24*f_{n-1} + 1/24*f{n-2})
% For the Cauchy problem dy/dx = -y - 3 * x with y(0) = 1, we need to solve
% a linear equation which is shown in the code
x(1) = xRange(1);
h = (xRange(2) - xRange(1)) / numSteps; %step size
y(:,1) = yInitial;

k = 1;
fValue1 = f_ode(x(k),y(:,k));
% Method1: Set the exact solution to the initial value
x(k+1) = x(k) + h;
y(:, k+1) = exmp_sol(x(k+1));
fValue2 = f_ode(x(k+1),y(:,k+1));

x(k+2) = x(k+1) + h;
y(:, k+2) = exmp_sol(x(k+2));

%Adam 3nd-order 
for k = 3:numSteps
    fValueold1 = fValue1;
    fValue1 = fValue2;
    fValueold2 = fValue2;
    fValue = f_ode(x(k),y(:,k));
    fValue2 = fValue;
    x(k+1) = x(k) + h;
    
    y(:,k+1) = y(:,k)/(1+3/8*h) + h/(1+3/8*h)*(-9/8*x(k+1)+19/24*fValue - 5/24*fValueold2 + 1/24*fValueold1);
end

