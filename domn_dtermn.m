%% Domain determination

function [t0, tend, tss] = domn_dtermn(t, x)
    figure(1); hold on;
    title("Raw Response Data");
    plot(t, x);
    xlabel("t");
    ylabel("x");

    t0 = input('Estimate the time the response starts: ');
    tend = input('Estimate the time stiction affects response: ');
    tss = input('Estimate the time response reaches steady-state: ');
    hold off;
    close 1;
    t0 = find(t >= t0);
    t0 = t0(1);
    tss = find(t <= tss);
    tss = tss(end);
    tend = find(t <= tend);
    tend = tend(end);
end
