%% iGEM Toronto 2017 Dry Lab MATLAB model %%
% Further documentation goes here.
%
%
%

%% Parameter Definitions %%
% In this section, the parameters for the model are listed and defined.
%
% time: this is the time vector for the model. It simply deflined by the
%       linspace function. The linspace funcation has the parameters (X1, X2, N),
%       where X1 is start point, X2 is end point, and N is the number of steps
%       between them.
time = linspace(0, 20, (20-0)/0.1);

% alpha_2:
alpha_2 = 0.92356;

% K_1:
K_1 = 0.3;

% R_1:
R_1 = 1.7;

% mu: this is the light function, defined as a square wave oscillating
%     between values 0 (off) and 1 (on). The duty cycle can be altered by
%     including factors to adjust the 't' value.
mu = @(t) (square(t, 50)+1)/2;

% psi_1: 
psi_1 = @(t) alpha_2/((K_1^n + R_1^n).^(1-mu(t)));

% gamma_2: 

gamma_2 = 1;

% alpha_theta:
alpha_theta = 1.6;

% gamma_theta:
gamma_theta = 1;

% k: 
k = alpha_theta/alpha_2;

% n: the exponent in the Michaelis Menten Equations
n = 1;

% alpha_x:
alpha_x = 1.4;

% gamma_L:
gamma_L = 0.6;

%% Listing of the Differential Equations %%
% In this section, the differential equations of interest for the model are
%   listed, defined, and sllved for. This model is currently solved using the
%   Runge-Kutta method, with the 'rk' function written by Peter Brinkmann. The
%   'rk2' function uses the same algorithm, but was slightly adapted to allow
%   for a third parameter to be included in the inline function.

% Equation 5 - dx2/dTau:
%   This is
dx2 = @(t, x2) psi_1(t) - gamma_2*x2;
solnX = rk(dx2, 1, time);

% Equation 6 - dTheta/dTau:
%   This is
dT = @(t, T) k.*psi_1(t) - gamma_theta*T;
solnT = rk(dT, 1, time);

% Equation 7 - dLambda/dTau:
%   This is
dL = @(t, L, x2) (alpha_x/(1 + x2^n)) - gamma_L*L;
solnL = rk2(dL, 0, solnX, time);

%% Plotting of the Solutions %%

plot(time, solnX, time, solnT, time, solnL, time, mu(time));
title('Plot of Equations 5, 6, and 7, along with Light Values');
ylabel('Concentrations');
xlabel('time');
legend('x2', 'theta', 'lambda', 'light', 'Location', 'NE');