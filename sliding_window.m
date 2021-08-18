function [f N] = sliding_window(x, y, xplot,w)
% computes sliding window
%
% [f N] = sliding_window(x, y, xplot, w)
%   x - x value of all datapoints
%   y - y value of all datapoints
%   xplot - x values to perform sliding window over
%   w - window size

for i=1:length(xplot)
    igood = find(x>xplot(i)-w & x<xplot(i)+w);
    if(length(igood)>3)
        f(i) = mean(y(igood));
    else
        f(i) = NaN;
    end
    N(i) = length(igood);
end
