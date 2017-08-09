function model_ODEs()
    % Solves the ODES given by the blue light paper
    % Reference: https://www.mathworks.com/help/symbolic/solve-a-system-of-differential-equations.html
    
    % represent differentials
    syms k_RFPmRNA(t) k_GFPmRNA(t) mRNA_RFP(t) mRNA_GFP(t) RFP(t) GFP(t);
    % represent parameters
    [k_sInd, L_I, d_sInd, k_sRep, d_sRep, L_R, k_basalmRNARFP, d_mRNARFP, ...
        k_basalmRNAGFP, d_mRNAGFP, k_RFP, k_GFP, d_RFP, d_GFP] = deal(0);
    
    % define the 6 equations
    ode1 = diff(k_RFPmRNA) == k_sInd*L_I - d_sInd*k_RFPmRNA;
    ode2 = diff(k_GFPmRNA) == k_sRep*L_R - d_sRep*k_GFPmRNA;
    ode3 = diff(mRNA_RFP) == k_RFPmRNA + k_basalmRNARFP - d_mRNARFP*mRNA_RFP;
    ode4 = diff(mRNA_GFP) == k_GFPmRNA + k_basalmRNAGFP - d_mRNAGFP*mRNA_GFP;
    ode5 = diff(RFP) == k_RFP*mRNA_RFP - d_RFP*RFP;
    ode6 = diff(GFP) == k_GFP*mRNA_GFP - d_GFP*GFP;
    odes = [ode1; ode2; ode3; ode4; ode5; ode6];

    % solve the ODEs with initial conditions (C1-C6)
    cond1 = k_RFPmRNA(0) == 0;
    cond2 = k_GFPmRNA(0) == 1;
    cond3 = mRNA_RFP(0) == 1;
    cond4 = mRNA_GFP(0) == 1;
    cond5 = RFP(0) == 0;
    cond6 = GFP(0) == 0;
    conds = [cond1; cond2; cond3; cond4; cond5; cond6];
    [k_RFPmRNA_Sol(t), k_GFPmRNA_Sol(t), mRNA_RFP_Sol(t), mRNA_GFP_Sol(t), ...
        RFP_Sol(t), GFP_Sol(t)] = dsolve(odes, conds);

    % Note: can simplify the equation by calling S = simplify(S)
    % but unnecessary in this case
    
    % plot
    fplot(k_RFPmRNA_Sol);
    hold on;
    fplot(k_GFPmRNA_Sol);
    hold on;
    fplot(mRNA_RFP_Sol);
    hold on;
    fplot(mRNA_GFP_Sol);
    hold on;
    fplot(RFP_Sol);
    hold on;
    fplot(GFP_Sol);
    grid on;
    legend('k_RFPmRNA_Sol','k_GFPmRNA_Sol', 'mRNA_RFP_Sol', 'mRNA_GFP_Sol', ...
        'RFP_Sol', 'GFP_Sol', 'Location','best');
    title ('Blue Light Mediated Transcriptional Activation');
    
    % plotting with different initial variables example:
    %for i = 0:10
    %    fplot(subs(ySol,'t',i),[0 1])
    %end
    %hold off
    %axis([0 1 -1 25])
    %title('Test Figure for dSolve')