%% Calculate system params
function [m, c, k] = param_solve(wn1, wn2, M1, M2, zeta)
	m = (M2*wn2^2 - M1*wn1^2)/(wn1^2-wn2^2);
    k = wn1^2*(m+M1);
    c = 2*zeta*sqrt(k*(m+M1));
end