% Ben Boggs - 2nd Order System Identification

clearvars; clc; close all;

%% Input: Define Step Magnitude

A = 1;

%% Data read in

d1 = csvread("free response.csv",1,0);
d2 = csvread("forced response.csv",1,0);

free = [d1(:,1),d1(:,10)];
forced = [d2(:,1),d2(:,10)];

%% Domain determination

figure(1); hold on;
title("Raw Free Response Data");
plot(free(:,1), free(:,2));
xlabel("t");
ylabel("x");

free_t0 = input('Estimate the time the response starts: ');
free_tend = input('Estimate the time stiction affects response: ');
free_ss = input('Estimate the time response reaches steady-state: ');
hold off;
close 1;

t = free(:,1);
t0 = find(t >= free_t0);
t0 = t0(1);
tss = find(t <= free_ss);
tss = tss(end);
free = [d1(t0:tss,1)-d1(t0,1), -d1(t0:tss,10)];
tend = find(t <= free_tend);
tend = tend(end);

figure(1); hold on;
title("Raw Free Response Data");
plot(forced(:,1), forced(:,2));
xlabel("t");
ylabel("x");

forced_t0 = input('Estimate the time the response starts: ');
forced_tend = input('Estimate the time stiction affects response: ');
forced_ss = input('Estimate the time response reaches steady-state: ');
hold off;
close 1;

t = forced(:,1);
t0 = find(t >= forced_t0);
t0 = t0(1);
tss = find(t <= forced_ss);
tss = tss(end);
forced = [d2(t0:tss,1)-d2(t0,1), d2(t0:tss,10)];
tfend = find(t <= forced_tend);
tfend = tfend(end);

%% Free Response Analysis

[~, locs] = findpeaks(free(1:tend,2));
free = [free(locs(1):end,1)-free(locs(1),1), free(locs(1):end,2)];
[pks, locs] = findpeaks(free(1:tend,2));
if locs(1) ~= 1
    pks = [free(1,2); pks];
    locs = [1;locs];
end
%figure(3); hold on;
%plot(free(:,1),free(:,2), free(locs,1), pks, '*');

figure(1); hold on;
plot(free(:,1),free(:,2));
title("Free Response");
xlabel("t");
ylabel("x");

[T, wd, dwd] = period(free(locs,1), length(locs));


% curve fit
curvefit = fit(free(locs,1), pks, 'exp1');
sigma = -curvefit.b;
x0 = curvefit.a;
%figure(10), hold on;
%plot(curvefit, free(:,1), free(:,2));

zeta = 1./sqrt(1+(2.*pi./sigma./T).^2);
dz = std(zeta);
zeta = mean(zeta)

wn = 2.*pi./(sqrt(1-zeta^2).*T);
dwn = std(wn);
wn = mean(wn)
t = free(:,1);
model = (x0/(1-zeta^2)).*exp(-sigma.*t).*cos(wd.*t);
pt = "Free Response: Curve Fit Method";
mplot(2, model, pt, t, free(:,2));

% logarithmic decrement

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
mplot(3, model, pt, free(:,1), free(:,2));

%% Forced Response Analysis

[pks, locs] = findpeaks(forced(1:tfend,2));
xss = forced(end,2)
%figure(3); hold on;
%plot(forced(:,1),forced(:,2), forced(locs,1), pks, '*');

figure(4); hold on;
plot(forced(:,1),forced(:,2));
plot(forced(:,1), xss.*ones(size(forced(:,1))), 'k--');
title("Step Response");
xlabel("t");
ylabel("x");

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
mplot(5, model, pt, t, forced(:,2));
plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');

% Rise Time
x = forced(:,2);
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
mplot(6, model, pt, t, x);
plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');

% curve fit
curvefit = fit(t(locs), pks-xss, 'exp1');
sigma = -curvefit.b;

%figure(10), hold on;
%plot(curvefit, t, x-xss);

zeta = 1./sqrt(1+(2.*pi./sigma./T).^2);
dz = std(zeta);
zeta = mean(zeta)
phi = atan2(sqrt(1-zeta^2),zeta)
wn = 2.*pi./(sqrt(1-zeta^2).*T);
dwn = std(wn);
wn = mean(wn)
model = xss.*(1-(exp(-sigma.*t).*sin(wd.*t+phi))./sqrt(1-zeta^2));
pt = "Step Response: Curve Fit Method";
mplot(7, model, pt, t, x);
plot(t, xss.*ones(size(t)), 'k--', 'HandleVisibility', 'off');

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