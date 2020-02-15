% Individual Sine Analysis

function [Mmax, halfpower, w, G, phi] = indsine(sinefolder)
    sinefiles = dir(fullfile(sinefolder, '*.csv'));
    for i = 1:length(sinefiles)
        baseFileName = sinefiles(i).name;
        filename = fullfile(sinefolder, baseFileName);
        nums = str2double(regexp(filename, '[\d\.]+', 'match'));
        w(i) = nums(1)*2*pi;
        d1 = csvread(filename,1,0);
        t = d1(:,1);
        x = d1(:,10);
        cx = d1(:,5);
        A = nums(2);
        [pks, locs] = findpeaks(x);
        G(i) = pks(end)/A;
        [~, clocs] = findpeaks(cx);
        [~,out,~,~,~] = filloutliers(x(locs),'previous');
        [~,outc,~,~,~] = filloutliers(cx(clocs),'previous');
        locs = locs(~out);
        clocs = clocs(~outc);
        s = t(clocs(end-5:end))-t(locs(end-5:end));
        phi(i) = w(i)*mean(s)/(2*pi)*360;
    end
    [w,I] = sort(w);
    G = G(I);
    phi = phi(I);
    wr = w(G == max(G));
    Mmax.wn = wr;
    halfpower.wn = wr;
    M = G/G(1);
    Mmax.zeta = 1/(2*max(M));
    hp = ones(size(w)).*(max(M)/sqrt(2));
    [wi,Mi] = polyxpoly(w,hp,w,M);
    halfpower.zeta = abs(wi(1)-wi(2))/(2*wr);
    
    wm = (0:1000);
    
    Mmax.model = (G(1)*wr^2)./sqrt((wr^2-wm.^2).^2 + (2*Mmax.zeta*wr.*wm).^2);
    halfpower.model = (G(1)*wr^2)./sqrt((wr^2-wm.^2).^2 + (2*halfpower.zeta*wr.*wm).^2);
    
    Mmax.phimodel = atan2((-2*Mmax.zeta*wr).*wm, (wr^2-wm.^2))./(2*pi).*360;
    halfpower.phimodel = atan2((-2*halfpower.zeta*wr).*wm, (wr^2-wm.^2))./(2*pi).*360;
end
