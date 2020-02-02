%% step response curve fit
function [wn, zeta, model] = curve_fit_step(t,x,tend,plt)
	forced = [t, x];
	[pks, locs] = findpeaks(forced(1:tend,2));
	xss = forced(end,2);
	curvefit = fit(t(locs), pks-xss, 'exp1');
	sigma = -curvefit.b;

	%figure(10), hold on;
	%plot(curvefit, t, x-xss);

    [T, wd, ~] = period(forced(locs,1), length(locs));
    zeta = 1./sqrt(1+(2.*pi./sigma./T).^2);
	dz = std(zeta);
	zeta = mean(zeta);
	phi = atan2(sqrt(1-zeta^2),zeta);
	wn = 2.*pi./(sqrt(1-zeta^2).*T);
	dwn = std(wn);
	wn = mean(wn);
	model = xss.*(1-(exp(-sigma.*t).*sin(wd.*t+phi))./sqrt(1-zeta^2));
	pt = "Step Response: Curve Fit Method";
	if plt > 0
		mplot(plt, model, pt, t, x, locs);
		plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');
	end
end