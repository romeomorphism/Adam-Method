% we plot the the error on the right hand value of the interval in log
% scale

n = 10;
num = [n]; % an array of number of grid points in the scheme
h = [2/n]; % an array of each stepsize for the array num
while n <= 1000
    n = 2*n;
    num(end+1) = n;
    h(end+1) = 2/n;
end

error1 = [];
error2 = [];
error3 = [];
error4 = [];
yInit = 1;

for numSteps = num
    [x1, y1] = Adam_explicit2( @exmp_ode, [0.0, 2.0], yInit, numSteps);
    [x2, y2] = Adam_explicit3( @exmp_ode, [0.0, 2.0], yInit, numSteps);
    [x3, y3] = Adam_implicit3( @exmp_ode, [0.0, 2.0], yInit, numSteps);
    [x4, y4] = Adam_implicit4( @exmp_ode, [0.0, 2.0], yInit, numSteps);
    error1(end+1) = abs(exmp_sol(2) - y1(end));
    error2(end+1) = abs(exmp_sol(2) - y2(end));
    error3(end+1) = abs(exmp_sol(2) - y3(end));
    error4(end+1) = abs(exmp_sol(2) - y4(end));
end

slope1 = polyfit(log(h), log(error1), 1);
slope2 = polyfit(log(h), log(error2), 1);
slope3 = polyfit(log(h), log(error3), 1);
slope4 = polyfit(log(h), log(error4), 1);

loglog(h, error1, 'DisplayName','Adam Explicit 2')
hold on
loglog(h, error2, 'DisplayName','Adam Explicit 3')
hold on 
loglog(h, error3, 'DisplayName','Adam Implicit 3')
hold on 
loglog(h, error4, 'DisplayName','Adam Implicit 4')
hold on
xlabel('stepsize h')
ylabel('error of solution at x=2')
hold off
legend
    


    