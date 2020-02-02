%% Data read in

function [t,x] = readData(filename)
    d1 = csvread(filename,1,0);
    t = d1(:,1);
    x = d1(:,10);
end
