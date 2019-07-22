function [a,b] = delta_to_interval(mu,delta)
% convert mu +/- delta% to interval [a,b]

a = mu - delta*abs(mu);
b = mu + delta*abs(mu);