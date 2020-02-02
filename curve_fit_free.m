% curve fit
function [x0, wn, zeta, model] = curve_fit_free(t,x,tend,plt)

	[~, locs] = findpeaks(x(1:tend));
	free = [t(locs(1):end)-t(locs(1)), x(locs(1):end)];
	[pks, locs] = findpeaks(free(1:tend,2));
	if locs(1) ~= 1
	    pks = [free(1,2); pks];
	    locs = [1;locs];
	end
	%figure(3); hold on;
	%plot(free(:,1),free(:,2), free(locs,1), pks, '*');

	%{
	figure(1); hold on;
	plot(free(:,1),free(:,2));
	title("Free Response");
	xlabel("t");
	ylabel("x");
	%}

	[T, wd, dwd] = period(free(locs,1), length(locs));


	% curve fit
	curvefit = fit(free(locs,1), pks, 'exp1');
	sigma = -curvefit.b;
	x0 = curvefit.a;
	%figure(10), hold on;
	%plot(curvefit, free(:,1), free(:,2));

	zeta = 1./sqrt(1+(2.*pi./sigma./T).^2);
	dz = std(zeta);
	zeta = mean(zeta);

	wn = 2.*pi./(sqrt(1-zeta^2).*T);
	dwn = std(wn);
	wn = mean(wn);
	t = free(:,1);
	model = (x0/(1-zeta^2)).*exp(-sigma.*t).*cos(wd.*t);
	pt = "Free Response: Curve Fit Method";
	if plt > 0
	    mplot(plt, model, pt, t, free(:,2), locs);
	end
end
