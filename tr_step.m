% Rise Time
function [wn, zeta, model] = tr_step(t,x,tend,plt)
	forced = [t, x];
	[pks, locs] = findpeaks(forced(1:tend,2));
	xss = forced(end,2);
	i2 = find(x >= xss);
	i2 = i2(1);
	i1 = i2-1;
	% interpolate (linear)
	tr = (xss-x(i1))*(t(i2)-t(i1))/(x(i2)-x(i1)) + t(i1);

    [~, wd, ~] = period(forced(locs,1), length(locs));
    
    phi = pi-tr*wd;
	zeta = sqrt(1/(tan(phi)^2+1));
	wn = wd/sqrt(1-zeta^2);
	sigma = zeta*wn;

	model = xss.*(1-(exp(-sigma.*t).*sin(wd.*t+phi))./sqrt(1-zeta^2));

	pt = "Step Response: Rise Time Method";
	if plt > 0
		mplot(plt, model, pt, t, x, locs);
		plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');
	end
end