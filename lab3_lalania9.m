%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Improved Euler 
% method (also known as the Heun method), and compare its results to those 
% of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in on the
% due date. Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online on Quercus.
%
% MAT292, Fall 2019, Stinchcombe & Parsch, modified from
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, modified from
% MAT292, Fall 2013, Sinnamon & Sousa, modified from
% MAT292, Fall 2011, Hart & Pym

%% Student Information
%
% Student Name: Armaan Lalani
%
% Student Number: 1005023225
%

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now MATLAB can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%

% NOTES: The function is at 'improvedEuler.m'.

%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the following initial-value problems (from lab 2, 
% exercises 1, 4-6), approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the graphs of your Improved Euler Approximation with the |ode45| 
% approximation.
%
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
%
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
%
% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
%
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
%
% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences for
% each of the four IVPs.

% ANSWER:
% Both parts a and b are very similar to the actual solution. For part c,
% there is a different between Euler's method and the actual solution in a
% certain region, but is relatively the same for the mejority of the
% solution. For part d, the solution diverges before t gets to 1, and the
% slution will stop around 0.5 because of a restriction placed by MATLAB.
% Euler's method will continue to compute until t=1 on the other hand.

clear; clc; close;

fa = @(t, y) y * tan(t) + sin(t);
ya0 = -0.5;
ta = [0, pi];

fb = @(t, y) 1/y^2;
yb0 = 1;
tb = [1, 10];

fc = @(t, y) 1 - t * y / 2;
yc0 = -1;
tc = [0, 10];

fd = @(t, y) y^3 - t^2;
yd0 = 1;
td = [0, 1];
figure;
hold on;

[t, y] = improvedEuler(fa, ta(1), ta(2), ya0, 0.001);
soln = ode45(fa, ta, ya0);

title('A');
plot(t, y, 'g-', 'LineWidth', 2);
plot(soln.x, soln.y, 'b--', 'LineWidth', 2);
legend('Improved Euler', 'ODE45', 'Location', 'Southeast');
xlabel('t');
ylabel('y(t)');
figure;
hold on;

[t, y] = improvedEuler(fb, tb(1), tb(2), yb0, 0.001);
soln = ode45(fb, tb, yb0);

title('B');
plot(t, y, 'g-', 'LineWidth', 2);
plot(soln.x, soln.y, 'b--', 'LineWidth', 2);
legend('Improved Euler', 'ODE45', 'Location', 'Southeast');
xlabel('t');
ylabel('y(t)');
figure;
hold on;

[t, y] = improvedEuler(fc, tc(1), tc(2), yc0, 0.001);
soln = ode45(fc, tc, yc0);

title('C');
plot(t, y, 'g-', 'LineWidth', 2);
plot(soln.x, soln.y, 'b--', 'LineWidth', 2);
legend('Improved Euler', 'ODE45', 'Location', 'Southeast');
xlabel('t');
ylabel('y(t)');
figure;
hold on;

[t, y] = improvedEuler(fd, td(1), td(2), yd0, 0.001);
soln = ode45(fd, td, yd0);

title('D');
plot(t, y, 'g-', 'LineWidth', 2);
plot(soln.x, soln.y, 'b--', 'LineWidth', 2);
legend('Improved Euler', 'ODE45', 'Location', 'Southeast');
xlabel('t');
ylabel('y(t)');

%% Exercise 3
%
% Objective: Use Euler's method and verify an estimate for the global error.
%
% Details: 
%
% (a) Use Euler's method (you can use
% euler.m from iode) to solve the IVP
%
% |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
%
% from |t=0| to |t=0.5|.
%
% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
%
% (c) Read the attached derivation of an estimate of the global error for 
%     Euler's method. Type out the resulting bound for En here in
%     a comment. Define each variable.
%
% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.
%
% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the order of Euler's method.

% ANSWERS:
% (b) 2tdt = dy/sqrt(1-y^2)
% t^2 + c = arcsin(y)
% y(t) = sin(t^2 + c)
% y(0) = 0
% y(t) = sin(t^2)

% y(0.5) = sin(0.25)

% (c) f_max = 2 * 0.5 * sqrt(1 - sin(0.5)^2) = 0.878
% ft_max = 2 * sqrt(1 - sin(0.5) ^ 2) = 1.755
% fy_max = 2 * 0.5 * sin(0.5) / sqrt(1 - sin(0.5)^2) = 0.546
% M ~ 1.755

% (d) for dt=0.01, dt*n=0.5 ==> En <= (1+M)*dt/2 * (exp(M*dt*n) - 1) ~ 0.0194
% 0.0047 was the actual error, which is less than the estimated error

% (e) For dt=0.01 the global error was ~0.0047, and for dt=0.1 the global
% error was ~0.048. As displayed, the error increases by a factor of 10
% when the step size increases by a factor of 10. This shows the
% correlation of error and dt is linear.

close all; clear; clc;

hold on;
title('Euler Method');

clc;
f = @(t, y) 2 * t * sqrt(1 - y^2);
y0 = 0;

step = 0.01;
c = 0 : step : 0.5;
y = euler(f, y0, c);
fprintf('Error with step=%g is %g\n', step, abs(sin(0.5^2) - y(end)));
plot(c, y, 'b-', 'LineWidth', 2);

step = 0.1;
c = 0 : step : 0.5;
y = euler(f, y0, c);
fprintf('Error with step=%g is %g\n', step, abs(sin(0.5^2) - y(end)));
plot(c, y, 'g-', 'LineWidth', 2);

plot(c, sin(c.^2), 'r--', 'LineWidth', 2);

xlabel('t');
ylabel('y(t)');
legend('0.01 step', '0.1 step', 'Actual Solution', 'Location', 'SouthEast');

%% Adaptive Step Size
%
% As mentioned in lab 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where |h| is an initial step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Create an implementation of Euler's method by modifying your solution to 
% exercise 1. Change it to include the following:
%
% (a) On each timestep, make two estimates of the value of the solution at
% the end of the timestep: |Y| from one Euler step of size |h| and |Z| 
% from two successive Euler steps of size |h/2|. The difference in these
% two values is an estimate for the error.
%
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be
% successful and set the new solution value to be |Z+D|. This value has
% local error |O(h^3)|. If |abs(D)>=tol|, reject this step and repeat it 
% with a new step size, from (c).
%
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
%
% Comment on what the formula for updating the step size is attempting to
% achieve.

% ANSWERS:
% Euler's method is more accurate based on the smaller step size that is
% used. The absolute value is used for indicating the error, so when the
% ratio is less than 1, the error is higher than the threshold value, which
% leads to the step size needing to be decreased to improve on this
% accuracy. When the ratio is greater than 1, the error value is acceptable
% and the function continues.

%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.
%
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
%
% (c) Plot both approximations together with the exact solution.

clc; clear; close all;

f = @(t, y) 2 * t * sqrt(1 - y^2);
t_euler = 0 : 0.025 : 0.75;
y_euler = euler(f, 0, t_euler);

[t_adaptive, y_adaptive] = adaptiveEuler(f, 0, 0.75, 0, 0.025);

hold on;
title('Adaptive Euler Method');
xlabel('t');
ylabel('y(t)');

plot(t_euler, y_euler, 'g-', 'LineWidth', 2);
plot(t_adaptive, y_adaptive, 'r-', 'LineWidth', 2);
plot(t_euler, sin(t_euler.^2), 'b--', 'LineWidth', 2);

legend('Euler', 'Adaptive Euler', 'Exact Solution', 'Location',...
    'SouthEast');

%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.
% 
% (b) Plot the exact solution (from exercise 3.b), the Euler's 
% approximation (from exercise 3.a) and the adaptive Euler's approximation 
% (from exercise 5) from |t=0| to |t=1.5|.
%
% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.

% ANSWERS:
% (a) Adaptive Euler's method is closer to the actual solution because this
% method adjusts the stepsize to accomodate for the error. If the error is
% above some threshold value, the stepsize is decreased in order to
% accomodate and improve on accuracy.

% (c) Regular Euler's method fails because the change in the values of f
% occur very quickly, so the error in this region are also very large.
% Adaptive Euler's Method is slightly better in this case, but the same
% thing occurs. The error indicator for adaptive euler ends up being much
% smaller than the actual error, which leads to the stepsize not being
% decreased enough to keep the euler solution close enough to the actual
% solution.

clc; clear; close all;

f = @(t, y) 2 * t * sqrt(1 - y^2);
t_euler = 0 : 0.025 : 1.5;
y_euler = euler(f, 0, t_euler);

[t_adaptive, y_adaptive] = adaptiveEuler(f, 0, 1.5, 0, 0.025);

hold on;
title('Adaptive Euler');
xlabel('t');
ylabel('y(t)');

plot(t_euler, y_euler, 'g-', 'LineWidth', 2);
plot(t_adaptive, y_adaptive, 'r-', 'LineWidth', 2);
plot(t_euler, sin(t_euler.^2), 'k--', 'LineWidth', 2);
legend('Euler', 'Adaptive Euler', 'Exact Solution', 'Location',...
    'Southeast');