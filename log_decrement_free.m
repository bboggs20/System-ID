% logarithmic decrement
function [wn, zeta, model] = log_decrement_free(t,x,tend,plt)

	[~, locs] = findpeaks(x(1:tend));
	free = [t(locs(1):end)-t(locs(1)), x(locs(1):end)];
    x = free(:,2);
    t = free(:,1);
    x0 = x(1);
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
		mplot(plt, model, pt, free(:,1), free(:,2), locs);
	end
end