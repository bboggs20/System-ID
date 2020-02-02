%% Period/wd determination

function [T, wd, dwd] = period(tpks, n_peaks)
    for i=1:n_peaks-1
        T(i) = tpks(i+1)-tpks(i);
    end
    wd = 2.*pi./T;
    dwd = std(wd);
    wd = mean(wd);
end