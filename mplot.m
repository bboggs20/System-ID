%% Model Overplot

function [] = mplot(fignum, model, pltitle, t, x, tpks)
    figure(fignum); hold on;
    title(pltitle);
    xlabel("t");
    ylabel("x");
    plot(t, x, t, model);
    plot(t(tpks), x(tpks), 'r*');
    legend("Experiment", "Model", "Peaks Analyzed");
end