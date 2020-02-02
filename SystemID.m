% Ben Boggs - 2nd Order System Identification

%% Find system parameters

function [sys_model] = SystemID(A, free, forced, M1, M2)
	% given cleaned data, find relevant parameters using all possible methods
    % params:
    %   A: step gain (V)
    %   free: struct with 6 fields: t1,x1,tend1 for M1 and t2,x2,tend2 for M2
    %   forced: same as free but with step response data instead of free response
    %   M1: value of M1 in g
    %   M2: value of M2 in g
    %
    % The output is a struct of structs containing system parameters. 
    % Additionally, the model for all M1 systems will be plotted over the data.
    
   

    % analyze free responses:

    %	log dec
    [x0,wn1,zeta,model] = log_decrement_free(free.t1, free.x1, free.tend1, 1);
    [~,wn2,~,~] = log_decrement_free(free.t2, free.x2, free.tend2, 0);
    [m, c, k] = param_solve(wn1, wn2, M1, M2, zeta);
    wn = wn1;
    sys_model.ld = struct('m', m, 'c', c, 'k', k, 'wn', wn, 'zeta', zeta, 'x0', x0, 'model', model);

    %	curve fit
    [x0,wn1,zeta,model] = curve_fit_free(free.t1, free.x1, free.tend1, 2);
    [~,wn2,~,~] = curve_fit_free(free.t2, free.x2, free.tend2, 0);
    [m, c, k] = param_solve(wn1, wn2, M1, M2, zeta);
    wn = wn1;
    sys_model.fcfit = struct('m', m, 'c', c, 'k', k, 'wn', wn, 'zeta', zeta, 'x0', x0, 'model', model);

    % analyze step responses:

    %	Max Overshoot
    [wn1,zeta,model] = Mp_step(forced.t1, forced.x1, forced.tend1, 3);
    [wn2,~,~] = Mp_step(forced.t2, forced.x2, forced.tend2, 0);
    [m, c, k] = param_solve(wn1, wn2, M1, M2, zeta);
    wn = wn1;
    K_hw_ol = forced.x1(end)*k/A;
    sys_model.Mp = struct('m', m, 'c', c, 'k', k, 'wn', wn, 'zeta', zeta, 'K_hw_ol', K_hw_ol, 'model', model);

    %	Rise Time
    [wn1,zeta,model] = tr_step(forced.t1, forced.x1, forced.tend1, 4);
    [wn2,~,~] = tr_step(forced.t2, forced.x2, forced.tend2, 0);
    [m, c, k] = param_solve(wn1, wn2, M1, M2, zeta);
    wn = wn1;
    K_hw_ol = forced.x1(end)*k/A;
    sys_model.tr = struct('m', m, 'c', c, 'k', k, 'wn', wn, 'zeta', zeta, 'K_hw_ol', K_hw_ol, 'model', model);

    %	Curve Fit
    [wn1,zeta,model] = curve_fit_step(forced.t1, forced.x1, forced.tend1, 5);
    [wn2,~,~] = curve_fit_step(forced.t2, forced.x2, forced.tend2, 0);
    [m, c, k] = param_solve(wn1, wn2, M1, M2, zeta);
    wn = wn1;
    K_hw_ol = forced.x1(end)*k/A;
    sys_model.scfit = struct('m', m, 'c', c, 'k', k, 'wn', wn, 'zeta', zeta, 'K_hw_ol', K_hw_ol, 'model', model);
end
