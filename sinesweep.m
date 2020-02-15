% Sine Sweep Analysis

function [Mmax, halfpower, w, G] = sinesweep(filename, A)
    [t,x] = readData(filename);

    [pks, locs] = findpeaks(x);
    periods = diff(t(locs));
    win = 2*pi./periods;
    
    f = fit(t(locs(2:end)), win, 'Poly1');
    m = f.p1; b = f.p2;

    % positive peaks only
    pks = pks(pks > 0);
    locs = locs(x(locs) > 0);
 
    G = pks/A;
    w = t(locs)*m + b;
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
end