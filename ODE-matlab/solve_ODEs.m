function solve_ODEs()
    % Solves the ODES given by the blue light paper
    % Reference: https://www.mathworks.com/help/symbolic/solve-a-system-of-differential-equations.html
    
    % represent variables
    syms u(t) v(t);
    
    % define equations
    ode1 = diff(u) == 3*u + 4*v;
    ode2 = diff(v) == -4*u + 3*v;
    odes = [ode1; ode2];
    
    % solve the equations
    cond1 = u(0) == 0; % set C1
    cond2 = v(0) == 1; % set C2
    conds = [cond1; cond2];
    [uSol(t), vSol(t)] = dsolve(odes,conds);

    % can simplify the equation by calling S = simplify(S)
    
    % plot
    fplot(uSol);
    hold on;
    fplot(vSol);
    grid on;
    legend('uSol','vSol','Location','best');
    title ('Test ODE Plot');
    
    % plotting with different initial variables example:
    %for i = 0:10
    %    fplot(subs(ySol,'t',i),[0 1])
    %end
    %hold off
    %axis([0 1 -1 25])
    %title('Test Figure for dSolve')