% Maximum Overshoot
function [wn, zeta, model] = Mp_step(t,x,tend,plt)
	forced = [t, x];
	[pks, locs] = findpeaks(forced(1:tend,2));
	xss = forced(end,2);
	%figure(3); hold on;
	%plot(forced(:,1),forced(:,2), forced(locs,1), pks, '*');
	%{
	figure(4); hold on;
	plot(forced(:,1),forced(:,2));
	plot(forced(:,1), xss.*ones(size(forced(:,1))), 'k--');
	title("Step Response");
	xlabel("t");
	ylabel("x");
	%}

	[T, wd, dwd] = period(forced(locs,1), length(locs));

	% % Maximum Overshoot
	Mp = (pks(1) - xss)/xss;
	zeta = abs(log(Mp)/sqrt(log(Mp)^2 + pi^2));
	wn = wd/sqrt(1-zeta^2);
	dwn = dwd/sqrt(1-zeta^2);
	sigma = zeta*wn;
	phi = atan2(sqrt(1-zeta^2),zeta);

	t = forced(:,1);
	model = xss.*(1-(exp(-sigma.*t).*sin(wd.*t+phi))./sqrt(1-zeta^2));

	pt = "Step Response: % Max Overshoot Method";
	if plt > 0
		mplot(plt, model, pt, t, forced(:,2), locs);
		plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');
	end
end
