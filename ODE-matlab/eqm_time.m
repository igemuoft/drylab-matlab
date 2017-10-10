function out = eqm_time(gamma_2, gamma_lambda, time, period);
    % Change K and Gamma_Theta
    N = length(gamma_lambda);
    out = zeros(N, N);

    % n: the exponent in the Michaelis Menten Equations
    n = 1;

    % alpha_theta:
    alpha_theta = 1.6;

    % alpha_2:
    alpha_2 = 0.92356;

    % alpha_x:
    alpha_x = 1.4;

    % k: 
    K = alpha_theta/alpha_2;

    % K_1:
    K_1 = 0.3;

    % R_1:
    R_1 = 1.7;

    % mu: this is the light function, defined as a square wave oscillating
    %     between values 0 (off) and 1 (on). The duty cycle can be altered by
    %     including factors to adjust the 't' value.
    %mu = @(t) (square(t, 50)+1)/2; % using square function (works)
    mu = @(t) mod(ceil(t./period),2); % w/o square() - also works now

    % psi_1: 
    psi_1 = @(t) alpha_2/((K_1^n + R_1^n).^(1-mu(t)));

    for i=1:N
        for j = 1:N
            % Equation 5 - dx2/dTau:
            dx2 = @(t, x2) psi_1(t) - gamma_2(i,j).*x2;
            X2 = rk(dx2, 1, time);

            % Equation 7 - dLambda/dTau:
            dL = @(t, L, x2) (alpha_x/(1 + x2^n)) - gamma_lambda(i,j).*L;
            L = rk2(dL, 0, X2, time);

            % psi_1 when light is on or off
            psi_1_on = alpha_2;
            psi_1_off = alpha_2/(K_1^n + R_1^n);
            
            % solutions for on:
            x2_on = psi_1_on ./ gamma_2(i,j);
            eqm_value = alpha_x / ((1+x2_on^n) * gamma_lambda(i,j));
            
            % solutions for off:
            x2_off = psi_1_off ./ gamma_2(i,j);
            init_value = alpha_x / ((1+x2_off^n) * gamma_lambda(i,j));
        
            for k = 1:N
                if((L(k)-init_value) >= (eqm_value-init_value)*0.99)
                    out(i,j) = time(k)*10;
                    break
                end
            end
            
        
        end


    end

end