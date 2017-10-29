function out = parasweep_G2_A2(gamma_2, alpha_2, ON);
% function out=parasweep_G2_A2(gamma_2, alpha_2, ON);
%
% parasweep_G2_A2: performs a parameter sweep by evaluating the solution to
% the differential equation for x_2 at differnt values of gamma_2 and
% alpha_2
%
% Usage: out=parasweep_G2_A2(gamma_2, alpha_2, 1);
%
% Parameters:
%	gamma_2: a matrix of values, ideally generated through the meshgrid function
%	alpha_2: a matrix of values, ideally generated through the meshgrid
%	function, the transpose of gamma_2
%   ON:  binary value for if the light is on or off
%	
% Returns:
%	out: matrix of solution values, each value associated with the
%	corresponding index of the input matrices.

N = length(gamma_2);
out = zeros(N, N);

start_time = 0;
end_time = 6;
precision = 0.1; % lower = more accurate
light_period = 2;

total_points = (end_time - start_time) / precision;
time = linspace(start_time, end_time, total_points);

% n: the exponent in the Michaelis Menten Equations
n = 1;

% % alpha_2:
% alpha_2 = 0.92356;

% K_1:
K_1 = 0.3;

% R_1:
R_1 = 1.7;

% mu: this is the light function, defined as a square wave oscillating
%     between values 0 (off) and 1 (on). The duty cycle can be altered by
%     including factors to adjust the 't' value.
%mu = @(t) (square(t, 50)+1)/2; % using square function (works)
mu = @(t) mod(ceil(t./light_period),2); % w/o square() - also works now

% plot(time, mu(time))

% psi_1: 
% psi_1 = @(t) alpha_2/((K_1^n + R_1^n).^(1-mu(t)));

    for i=1:N

        for j = 1:N
            
            
            psi_1_on = alpha_2(i,j);
            psi_1_off = alpha_2(i,j)./(K_1^n + R_1^n);

            x2_on = psi_1_on ./ gamma_2(i,j);
            x2_off = psi_1_off ./ gamma_2(i,j);

            psi_1 = @(t) alpha_2(i,j)./((K_1^n + R_1^n).^(1-mu(t)));
            % Equation 5 - dx2/dTau:
            dx2 = @(t, x2) psi_1(t) - gamma_2(i,j).*x2;
            solnX = rk(dx2, 1, time);

            if(ON)
                out(i, j) = solnX(end_time/precision);
            else
                out(i, j) = solnX((2*light_period-1)/precision);
            end
        end


    end

end