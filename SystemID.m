% Ben Boggs - 2nd Order System Identification


function [] = SystemID()
end

%% Data read in

function [t,x] = readData(filename)
    d1 = csvread(filename,1,0);
    t = d1(:,1);
    x = d1(:,10);
end
%% Domain determination

function [t0, tend, tss, tp, xp] = domn_dtermn(t, x)
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
    tp = t(t0:tss)-t(t0);
    xp = x(t0:tss);
end

%% Free Response Analysis

% curve fit
function [T, wd, wn, zeta, model] = curve_fit_free(t,x,tend, plt)

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
	    mplot(plt, model, pt, t, free(:,2));
	end
end

% logarithmic decrement
function [T, wd, wn, zeta, model] = log_decrement_free(t,x,tend,plt)

	[~, locs] = findpeaks(x(1:tend));
	free = [t(locs(1):end)-t(locs(1)), x(locs(1):end)];
	[pks, locs] = findpeaks(free(1:tend,2));
	if locs(1) ~= 1
	    pks = [free(1,2); pks];
	    locs = [1;locs];
	end

	[T, wd, dwd] = period(free(locs,1), length(locs));

	for i=2:length(locs)
	    delta(i-1) = log(pks(1)/pks(i))/i;
	end
	zeta = mean(delta)/(2*pi);
	dz = std(delta./(2*pi));

	wn = 2.*pi./(sqrt(1-zeta^2).*T);
	dwn = std(wn);
	wn = mean(wn);

	sigma = zeta*wn;

	model = (x0/(1-zeta^2)).*exp(-sigma.*t).*cos(wd.*t);
	pt = "Free Response: Logarithmic Decrement Method";
	if plt > 0
		mplot(plt, model, pt, free(:,1), free(:,2));
	end
end

%% Forced Response Analysis

% Maximum Overshoot
function [T, wd, wn, zeta, model] = Mp_step(t,x,tend,plt)
	forced = [t, x];
	[pks, locs] = findpeaks(forced(1:tend,2));
	xss = forced(end,2)
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
		mplot(plt, model, pt, t, forced(:,2));
		plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');
	end
end

% Rise Time
function [T, wd, wn, zeta, model] = tr_step(t,x,tend,plt)
	forced = [t, x];
	[pks, locs] = findpeaks(forced(1:tend,2));
	xss = forced(end,2);
	i2 = find(x >= xss);
	i2 = i2(1);
	i1 = i2-1;
	% interpolate (linear)
	tr = (xss-x(i1))*(t(i2)-t(i1))/(x(i2)-x(i1)) + t(i1);

	phi = pi-tr*wd;
	zeta = sqrt(1/(tan(phi)^2+1));
	wn = wd/sqrt(1-zeta^2);
	sigma = zeta*wn;

	model = xss.*(1-(exp(-sigma.*t).*sin(wd.*t+phi))./sqrt(1-zeta^2));

	pt = "Step Response: Rise Time Method";
	if plt > 0
		mplot(plt, model, pt, t, x);
		plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');
	end
end

% curve fit
function [T, wd, wn, zeta, model] = curve_fit_step(t,x,tend,plt)
	forced = [t, x];
	[pks, locs] = findpeaks(forced(1:tend,2));
	xss = forced(end,2);
	curvefit = fit(t(locs), pks-xss, 'exp1');
	sigma = -curvefit.b;

	%figure(10), hold on;
	%plot(curvefit, t, x-xss);

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
		mplot(plt, model, pt, t, x);
		plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');
	end
end

%% Period/wd determination

function [T, wd, dwd] = period(tpks, n_peaks)
    for i=1:n_peaks-1
        T(i) = tpks(i+1)-tpks(i);
    end
    wd = 2.*pi./T;
    dwd = std(wd);
    wd = mean(wd);
end

%% Model Overplot

function [] = mplot(fignum, model, pltitle, t, x)
    figure(fignum); hold on;
    title(pltitle);
    xlabel("t");
    ylabel("x");
    plot(t, x, t, model);
    legend("Experiment", "Model");
end