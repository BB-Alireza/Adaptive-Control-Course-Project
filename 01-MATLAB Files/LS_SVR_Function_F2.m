%% Application of SUPPORT VECTOR REGRESSION to estimate an unknown non linear function of EV Model
% Least Squares Support Vector Regression method
clc;
clear all;
close all;

%% defining system parameters
Ra_plus_Rf = 0.12;      % Armature and Field Resistance
La_plus_Lf = 6.008e-3;  % Armature and Field Inductance
r = 0.25;               % Wheel Radious
J = 0.05;               % Inertai of the Motor
L_af = 1.766e-3;        % Mutual Inductance
I = 78;                 % Maximum Current
omega_norm = 2800;      % Nominal omega
m = 800;                % Total mass (kg)
A = 1.8;                % Frontal Area of the Vegicle
ro = 1.25;              % Air density
Cd = 0.3;               % Drag Coefficient
phi = 0;                % Hill climbing Angle
mu_rr = 0.015;          % Rolling Resistance Coefficient
G = 11;                 % Gearbox Ratio
B = 0.0002;             % Viscous Coefficient
g = 9.8;                % Gravity Constant 
V = 48;                 % Input Voltage  0 ~ 48 volt
%% Input and Output definition
%%% Initial Condition
x1 = 3;
x2 = 1.2;
x1_dot = (1/La_plus_Lf) * (-(Ra_plus_Rf) * x1 - L_af * x1 * x2) + (1/La_plus_Lf) * V;
x2_dot = (1 / (J + m * (r/G)^2)) * (L_af * x1^2 -B * x2 - (r/G) * (mu_rr * m * g + 0.5 * ro * A * Cd * (r/G)^2 * x2^2 + m*g*sin(phi)));

% Define the input range x.
x = [1:150]'; 
x_lenght=length(x);

%Define the original function to be estimated
y_orig= - 2*(Ra_plus_Rf) - 0.5*sin(2 *Ra_plus_Rf*x)+0.8*cos(4*x)-0.05*exp(x/30);

%Define the original function affected by noise that will generate the
%samples
y_noise=y_orig+0.3*randn(x_lenght,1);

% Define the training data for x
xs = [1:0.1:150]';
xs_lenght=length(xs);

% Define the training data for y
y_train=- 2*(Ra_plus_Rf) -0.5*sin(2 *Ra_plus_Rf*xs)+0.8*cos(4*xs)-0.05*exp(xs/30)+0.3*randn(xs_lenght,1);

% Create an empty vector to hold the approximate function values.
ys = zeros(size(x));

%% Start the Learning Algorithm, LEAST SQUARES-SVR

A = zeros(xs_lenght,xs_lenght); % initialize matrix A
C = 100; %Parameter defined to avoid overfitting
g = 0.01;%Radial Basis Functon learning parameter, is equal to 1/2sigma^2

%Radial Basis function
for i=1:xs_lenght
    for j=1:xs_lenght
        A(i,j)=exp(-g*(xs(i,1)-xs(j,1))^2);
        if i==j
            A(i,j)= A(i,j)+1/C;
        end
    end
end

O=[0,                 ones(1,xs_lenght);
   ones(xs_lenght,1), A];

b=zeros(xs_lenght+1,1);
c=zeros(xs_lenght+1,1);

 for f=1:xs_lenght
     c(f+1,1)=y_train(f,1);
 end 

b=inv(O)*c; %it contains LaGrange multipliers and bias

% For convenience, we separate the LaGrange Multipliers and the Bias term

%Bias term
bias=b(1,1);

%Define LaGrange Multipliers alfa
alfa=zeros(x_lenght,1); 

for i=1:xs_lenght
    alfa(i,1)=b(i+1,1); 
end

q=zeros(xs_lenght,1);% Storage variable

for j=1:x_lenght
   
 for i=1:xs_lenght
    
     q(i,1) = alfa(i,1).*exp(-g*(x(j,1)-xs(i,1))^2);
   
 end
     ys(j,1) = sum(q) + bias;
end

figure(1);
hold on; 

% Plot the original function as a blue line.
plot(x, y_orig,'color','blue',LineWidth=1.5);

% Plot the noisy data as red poinyts
plot(x, y_noise, '--','color','black',LineWidth=1.1);

% Plot the approximated function as a red line.
plot(x, ys,'color','red',LineWidth=1.5);

grid on
xlabel('Number of training samples')
legend('Original f_{2} Function', 'Noisy Samples', 'Approximated f_{2} Function');
title('LS-SVR Regression');

