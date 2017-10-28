function out = parasweep_kG(k, gamma_theta, ON);
% function out=parasweep_kG(k, gamma_theta, ON);
%
% parasweep_kG: performs a parameter sweep by evaluating the solution to
% the differential equation for 'theta' at differnt values of k and gamma_theta
%
% Usage: out=parasweep_kG(k, gamma_theta, 1);
%
% Parameters:
%	k: a matrix of values, ideally generated through the meshgrid function
%	gamma_theta: a matrix of values, ideally generated through the meshgrid
%	function, the transpose of k
%   ON:  binary value for if the light is on or off
%	
% Returns:
%	out: matrix of solution values, each value associated with the
%	corresponding index of the input matrices.

N = length(k);
out = zeros(N, N);

start_time = 0;
end_time = 15;
precision = 0.1; % lower = more accurate
light_period = 5;

total_points = (end_time - start_time) / precision;
time = linspace(start_time, end_time, total_points);

% n: the exponent in the Michaelis Menten Equations
n = 1;

% alpha_2:
alpha_2 = 0.92356;

% K_1:
K_1 = 0.3;

% R_1:
R_1 = 1.7;

% mu: this is the light function, defined as a square wave oscillating
%     between values 0 (off) and 1 (on). The duty cycle can be altered by
%     including factors to adjust the 't' value.
%mu = @(t) (square(t, 50)+1)/2; % using square function (works)
mu = @(t) mod(ceil(t./light_period),2); % w/o square() - also works now

% psi_1: 
psi_1 = @(t) alpha_2/((K_1^n + R_1^n).^(1-mu(t)));

    for i=1:N

        for j = 1:N


        % Equation 6 - dTheta/dTau:
        dT = @(t, T) k(i,j).*psi_1(t) - gamma_theta(i,j)*T;
        solnT = rk(dT, 1, time);

            if(ON)
                out(i, j) = solnT(end_time/precision);
            else
                out(i, j) = solnT((2*light_period-1)/precision);
            end
        end


    end

end