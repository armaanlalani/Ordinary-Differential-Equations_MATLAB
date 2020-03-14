%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
%
% This lab will teach you to solve ODEs using a built in MATLAB Laplace 
% transform function |laplace|. Also in this lab, you will write your own
% ODE solver using Laplace transforms and check whether the result yields
% the correct answer.
%
% You will learn how to use the |laplace| routine. 
% 
% There are five (5) exercises in this lab that are to be handed in.  
% Write your solutions in the template, including appropriate descriptions 
% in each step. Save the m-file and submit it on Quercus.
%
% Include your name and student number in the submitted file.
%
% MAT292, Fall 2019, Stinchcombe & Parsch, modified from
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, based on 
% MAT292, Fall 2013, Sinnamon & Sousa

%% Student Information
%
%  Student Name: Armaan Lalani
%
%  Student Number: 1005023225
%

%% Using symbolic variables to define functions
% 
% Recall the use of symbolic variables and function explained in the MATLAB
% assignment #2.

syms t s x y

f = cos(t)
h = exp(2*x)


%% Laplace transform and its inverse

% The routine |laplace| computes the Laplace transform of a function

F=laplace(f)

%%
% By default it uses the variable |s| for the Laplace transform
% But we can specify which variable we want:

H=laplace(h)
laplace(h,y)

% Observe that the results are identical: one in the variable |s| and the
% other in the variable |y|

%% 
% We can also specify which variable to use to compute the Laplace
% transform:

j = exp(x*t)
laplace(j)
laplace(j,x,s)

% By default, MATLAB assumes that the Laplace transform is to be computed
% using the variable |t|, unless we specify that we should use the variable
% |x|

%% 
% We can also use inline functions with |laplace|. When using inline
% functions, we always have to specify the variable of the function.

l = @(t) t^2+t+1
laplace(l(t))

%% 
% MATLAB also has the routine |ilaplace| to compute the inverse Laplace
% transform

ilaplace(F)
ilaplace(H)
ilaplace(laplace(f))

%% 
% If |laplace| cannot compute the Laplace transform, it returns an
% unevaluated call.

g = 1/sqrt(t^2+1)
G = laplace(g)

%% 
% But MATLAB "knows" that it is supposed to be a Laplace transform of a
% function. So if we compute the inverse Laplace transform, we obtain the
% original function

ilaplace(G)

%%
% The Laplace transform of a function is related to the Laplace transform 
% of its derivative:

syms g(t)
laplace(diff(g,t),t,s)


%% Exercise 1
%
% Objective: Compute the Laplace transform and use it to show that MATLAB
% 'knows' some of its properties.
%
% Details:  
%
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace
%   transform |F(s)|.
% (b) Find a function |f(t)| such that its Laplace transform is
%   |(s - 1)*(s - 2))/(s*(s + 2)*(s - 3)|
% (c) Show that MATLAB 'knows' that if |F(s)| is the Laplace transform of
%   |f(t)|, then the Laplace transform of |exp(at)f(t)| is |F(s-a)| 
% 
% (in your answer, explain part (c) using comments).      
%
% Observe that MATLAB splits the rational function automatically when
% solving the inverse Laplace transform.

clear all
syms t s x y

% (a)
f = exp(2*t)*t^3;
F = laplace(f); 

% F = 6/(s - 2)^4

% (b)
G = (s - 1)*(s - 2)/(s*(s + 2)*(s - 3));
g = ilaplace(G);

% g = (6*exp(-2*t))/5 + (2*exp(3*t))/15 - 1/3

% (c) 

syms f(t) a                 % F(s) is the Laplace transform of f(t)
F = laplace(f)              % F = laplace(f(t), t, s) = F(s)
G = laplace(exp(a*t)*f)     % G = laplace(f(t), t, s - a) = F(s-a)

%% Heaviside and Dirac functions
%
% These two functions are builtin to MATLAB: |heaviside| is the Heaviside
% function |u_0(t)| at |0|
%
% To define |u_2(t)|, we need to write

f=heaviside(t-2)
ezplot(f,[-1,5])

% The Dirac delta function (at |0|) is also defined with the routine |dirac|

g = dirac(t-3)

% MATLAB "knows" how to compute the Laplace transform of these functions

laplace(f)
laplace(g)


%% Exercise 2
%
% Objective: Find a formula comparing the Laplace transform of a 
%   translation of |f(t)| by |t-a| with the Laplace transform of |f(t)|
%
% Details:  
%
% * Give a value to |a|
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| 
%   and |F(s)| is the Laplace transform of |f(t)|, then find a 
%   formula relating |G(s)| and |F(s)|
%
% In your answer, explain the 'proof' using comments.

a = 2;
syms f(t)
g = heaviside(t-a)*f(t-a);
F = laplace(f)              % F = laplace(f(t), t, s)
G = laplace(g)              % G = exp(-2*s)*laplace(f(t), t, s)

% g(t) is a translation of f(t) by the value of a, and laplace(g) is
% laplace(f) multiplied by exp(-2*s) as shown above in the comments.

%% Solving IVPs using Laplace transforms
%
% Consider the following IVP, |y''-3y = 5t| with the initial
% conditions |y(0)=1| and |y'(0)=2|.
% We can use MATLAB to solve this problem using Laplace transforms:

% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y s

% Then we define the ODE

ODE=diff(y(t),t,2)-3*y(t)-5*t == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),1)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),2)

% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])

% We can check that this is indeed the solution

diff(y,t,2)-3*y


%% Exercise 3
%
% Objective: Solve an IVP using the Laplace transform
%
% Details: Explain your steps using comments
%
%
% * Solve the IVP
% *   |y'''+2y''+y'+2*y=-cos(t)|
% *   |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
% * for |t| in |[0,10*pi]|
% * Is there an initial condition for which |y| remains bounded as |t| goes to infinity? If so, find it.

% Define the unknown function and its corresponding variables
syms y(t) t Y s

% Right side of the function and its associated laplace transform
f = -cos(t);
F = laplace(f,t,s);

% Laplace transform of y'(t): Y1 = sY - y(0)
Y1 = s*Y;

% Laplace transform of y''(t): Y2 = sY1 - y'(0)
Y2 = s*Y1;

% Laplace transform of y'''(t): Y3 = sY2 - y''(0)
Y3 = s*Y2;

% Laplace transform of the left hand side equals to the right hand side set
% to 0, and solve for Y
Y = solve(Y3 + 2*Y2 + Y1 + 2*Y - F, Y);

% Inverse Laplace transform of the solution
y = ilaplace(Y,s,t);

ezplot(y,[0,300]);

%% Exercise 4
%
% Objective: Solve an IVP using the Laplace transform
%
% Details:  
% 
% * Define 
% *   |g(t) = 3 if 0 < t < 2|
% *   |g(t) = t+1 if 2 < t < 5|
% *   |g(t) = 5 if t > 5|
%
% * Solve the IVP
% *   |y''+2y'+5y=g(t)|
% *   |y(0)=2 and y'(0)=1|
%
% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.
%
% In your answer, explain your steps using comments.

% define the required variables
syms s t Y

% Write the right hand side using the heavyside function
f = 3*heaviside(t) + (t-2)*heaviside(t-2) + (4-t)*heaviside(t-5);

% Laplace transform of the right hand side
F =  laplace(f,t,s);

% Laplace transform of y'(t) : Y1 = s Y - y(0)
Y1 = s*Y - 2;

% Laplace transform of y''(t) : Y2 = s Y1 - y'(0)
Y2 = s*Y1 - 1;

% Set LHS - RHS = 0 and solve for Y
Y = solve(Y2 + 2*Y1 + 5*Y - F, Y);

% Inverse laplace transform of Y
y = ilaplace(Y, s, t);

ezplot(y,[0,12,0,2.25]);
%% Exercise 5a
%
% Objective: Use the Laplace transform to solve an integral equation
% 
% Verify that MATLAB knowns about the convolution theorem by explaining why the following transform is computed correctly.

syms t tau y(tau) s
I=int(exp(-2*(t-tau))*y(tau),tau,0,t)
laplace(I,t,s)

% The laplace transform of f*y=int(f(t-tau)y(tau)dtau) from 0 to t is F(s)Y(s). f(t) in this case is
% e^(-2t) which means that F(s) = 1/(s+2). Matlab returned the result laplace(y)/(s+2) which is
% the expected result.

%% Exercise 5b
% A particular machine in a factory fails randomly and needs to be replaced. Suppose that the times |t>=0| between failures are independent and identically distributed with probability density function |f(t)|. The mean number of failures |m(t)| at time |t| satisfies the renewal equation |m(t) = \int_0^t [1+m(t-tau)] f(tau) dtau|
%
% Details:  
%
% * Explain why the mean number of failures satisfies this intergal equation. Note that |m(0)=0|.
% * Solve the renewal equation for |m(t)| using MATLAB symbolic computation in the cases of i) exponential failure times |f(t) = exp(-t)| and ii) gamma-distributed failure times |f(t) = t^k/(k-1)! exp(-t)| for natural number |k|. Why does MATLAB have difficulty with the calculation for |k>=5|?
% * Verify the elementary renewal theorem: |m(t)/t| approaches the reciprocal of the mean of |f(t)| as |t| goes to infinity. 

% m(0) = 0 because no failure occurs before t=0.
% If we assume the value for m(t) is known for all 0<=t<T, then m(T) can be
% calculated. f(tau) is just the probability density function
% and m(T) are the expected values of failures.
% m(T) = \int_0^T x * f(tau).dtau
% where x is the number of failures occuring from tau and T,
% The expected number of failures occuring between tau and T is m(T-tau).
% Since there is a
% failure at t=tau, the total number of failures is m(T-tau)+1, which
% means x is just m(T-tau)+1.

syms f(t) m1(t) m2(t) t tau M1 M2

% Exponential decay probability density
f = exp(-t);

% Entering the equation
e = m1(t) - int((m1(t - tau) + 1) * subs(f, t, tau), tau, 0, t) == 0;

% Laplace transform
l = laplace(e);

% Substituting M and solving for M
l = subs(l, laplace(m1), M1);
M1 = solve(l, M1);

% inverse laplace transform. m=t
m1 = ilaplace(M1)

% Calculating the average of f(t) as t approaches infinity
average = eval(int(t * f, t, 0, inf));

% Multiplying m1(t) / t by average at t=infinity
result = subs(m1 / t, t, inf) * average;
result = eval(result)

% m(t) = t, so m(t)/t = 1. f(infinity) = e^(-inf) = 0.

% Gamma Distribution for the special case k=5
k = 5;
f = t^(k-1)/factorial(k-1) * exp(-t);

% New distribution
e = m2(t) == int((m2(t - tau) + 1) * subs(f, t, tau), tau, 0, t);
l = simplify(laplace(e));
l = simplify(subs(l, laplace(m2), M2));

M2 = simplify(solve(l, M2));
m2 = vpa(ilaplace(M2))

% Verifying the elementary renewal theorem
% Calculating the average of f(t) as t->infinity.
average = eval(int(t * f, t, 0, inf));

% Multiplying m(t) / t by average at t=infinity. 
result = subs(m2 / t, t, 1e10) * average;
result = eval(result)

% The reason that there is difficulty for k>=5 the process is slowed down:
% By solving this problem, M=k / (s(s+1)^k - ks). By using partial
% fractions, solving becomes difficult as the denominator's order increases.
% Hence, as k increases calculating ilaplace(M) gets more timely costly.
