%=============================================================================
% Project : Solving ODE With Different Numerical Methods
% Author  : Ritesh Tiwari
% Course  : MATLAB Online Training (Internshala) - Final Project
%
% Equation : (dy/dx) = (1 + 4*x) * sqrt(y) , Initial Condition: [x = 0, y = 1]
%
% Methods Implemented: 
%   1. Euler's Method
%   2. Heun's Method
%   3. Runge - Kutta 4th order Method
%   4. MATLAB Built-in Function[ode45()]
%   5. Analytical Method
%=============================================================================

clc, clearvars, close all, format compact

% Creating 'n' even spacing of x_values from [0 to 2]
n = input("enter the number of intervals[eg. 5, 10, 20]:");    % Default[n = 5]
x = linspace(0, 2, n);

% Creating input for choice between [Euler vs. Analytical, Heun vs. Analytical, Runge-Kutta vs Analytical, All (Euler + Heun + Runge-Kutta + ode45() + Analytical)]

fprintf("\nChoose a method to compare with analytical solution:\n")

fprintf("1 for Euler vs. Analytical\n2 for Heun vs. Analytical\n3 for Runge-Kutta vs Analytical\n4 for All (Euler + Heun + Runge-Kutta + ode45() + Analytical)\n\n")

choice = input("Enter your choice:");

% Defining a Variable 'x_an' for a Smoother Curve
x_an = linspace(0,2,41);    % Taken[n = 50]


% Custom Anonymous Function of the Respective given ODE 
f = @(x, y) (1 + 4*x) * sqrt(y);


% defining array for y_values of size 'n' for Euler's Method
E = linspace(0,0,n);

% defining array for y_values of size 'n' for Huen's Method
H = linspace(0,0,n);

% defining array for y_values of size 'n' for Runga-Kutta Method
RK = linspace(0,0,n);

% defining array for y_values of size 'n' for Analytical Method
An = linspace(0,0,(length(x_an)-1));


% Initial Value for all Methods
E(1) = 1; H(1) = 1; RK(1) = 1; An(1) = 1;

%=================================
% Analytical Method
%=================================
for i = 1:(length(x_an)-1)
    An(i+1) = (x_an(i+1)^2 + (1/2)*x_an(i+1) + 1).^2; % Main Equation
end

%=================================
% Euler's Method
%=================================
for i = 1:n-1
    h = (x(i+1) - x(i));                % Step Size - difference between two 'x' values[Eg: x(2) - x(1)]
    E(i+1) = E(i) + h * f(x(i), E(i));  % Main Equation
end

%=================================
% Huen's Method
%=================================
for i = 1:n-1
    S_left = f(x(i),H(i));                          % Left Slope
    h = (x(i+1) - x(i));                            % Step Size
    A = x(i) + h;                                   
    B = H(i) + h * S_left;
    S_right = f(A,B);                               % Right Slope
    H(i+1) = H(i) + (1/2) * h * (S_left + S_right); % Main Equation
end

%=================================
% Runga-Kutta Method
%=================================
for i = 1:n-1
    h = (x(i+1) - x(i));                                % Step Size
    k1 = h * f(x(i), RK(i));
    k2 = h * f(x(i) + (1/2)*h, RK(i) + (1/2)*k1);
    k3 = h * f(x(i) + (1/2)*h, RK(i) + (1/2)*k2);
    k4 = h * f(x(i) + h, RK(i) + k3);
    RK(i+1) = RK(i) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);  % Main Equation
end

%=================================
% MATLAB Built-in Method
%=================================
[t,y] = ode45(f,[0 2],1);   % Main Equation 

%=================================
% Plot Results
%=================================
figure;
hold on; grid on;

if choice == 1          % Euler vs. Analytical Method
    disp("Graph Produced Successfully")
    plot(x_an,An,"Color",'k','LineStyle','-','LineWidth',1.5)
    hold on
    plot(x,E,"Color",'r','Marker','v','LineStyle','--')
    legend('Analytical','Euler','Location','northwest'), xlabel('x'), ylabel('y')
    title('Solve ODE with Different Methods')

elseif choice == 2      % Huen vs. Analytical Method
    disp("Graph Produced Successfully")
    plot(x_an,An,"Color",'k','LineStyle','-','LineWidth',1.5)
    hold on
    plot(x,H,"Color",'g','Marker','o','LineStyle','--')
    legend('Analytical','Heun','Location','northwest'), xlabel('x'), ylabel('y')
    title('Solve ODE with Different Methods')

elseif choice == 3      % Runga-Kutta vs. Analytical Method
    disp("Graph Produced Successfully")
    plot(x_an,An,"Color",'k','LineStyle','-','LineWidth',1.5)
    hold on
    plot(x,RK,"Color",'b','Marker','square','LineStyle','--')
    legend('Analytical','Range-Kutta','Location','northwest'), xlabel('x'), ylabel('y')
    title('Solve ODE with Different Methods')

elseif choice == 4      % All Methods
    disp("Graph Produced Successfully")
    plot(x_an,An,"Color",'k','LineStyle','-','LineWidth',1.5)
    hold on
    plot(x,E,"Color",'r','Marker','v','LineStyle','--')
    hold on
    plot(x,H,"Color",'g','Marker','o','LineStyle','--')
    hold on
    plot(x,RK,"Color",'b','Marker','square','LineStyle','--')
    hold on
    plot(t',y',"Color",'m','Marker','pentagram','LineStyle','--')
    legend('Analytical','Euler','Heun','Range-Kutta','ode45','Location','northwest'), xlabel('x'), ylabel('y')
    title('Solve ODE with Different Methods')
 
else
    disp("INVALID INPUT --> CANNOT PRODUCE GRAPH!!!")

end