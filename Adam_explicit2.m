function[x, y] = Adam_explicit2(f_ode, xRange, yInitial, numSteps)
% Adam second-order method to solve ode dy/dx = f_ode(x, y)
% xRange = [x1, x2] where the solution is approximated in interval [x1, x2]
% yInitial = column vector of initial values for y at x1
% numSteps = number of equal stepsizes from x1 to x2
% x = row vector of value of x whose k-th elments is x1 + (k-1)*stepsize
% y = matrix whose k-th column is the approximate solution at x(k)

x(1) = xRange(1);
h = (xRange(2) - xRange(1)) / numSteps; %step size
y(:,1) = yInitial;

k = 1;
fValue = f_ode(x(k),y(k));
% Method1: Set the exact solution to the initial value
x(k+1) = x(k) + h;
y(:, k+1) = exmp_sol(x(k+1));

% Method2: Use onestep method to provide approximate solution to the initial value of multistep method
% Use RK2(Euler Halfstep) to approximate the solution at x(2)

%xhalf = x(k) + 0.5 * h;
%yhalf = y(:,k) + 0.5 * h * fValue;
%fValuehalf = f_ode(xhalf, yhalf);

%x(k+1) = x(k) + h;
%y(:,k+1) = y(:,k) + h * fValuehalf;

%Adam 2-order 
for k = 2:numSteps
    fValueold = fValue;
    fValue = f_ode(x(k), y(:,k));
    x(k+1) = x(k) + h;
    y(:,k+1) = y(:,k) + h * (3 * fValue - fValueold) / 2;
end





