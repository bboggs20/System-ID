%% data cut

function [tn, xn] = clean(t, x, t0, tss)
    tn = t(t0:tss)-t(t0);
    xn = x(t0:tss);
end
