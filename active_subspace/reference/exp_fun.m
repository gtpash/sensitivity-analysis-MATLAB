function y = exp_fun(x)

if isstruct(x)==1
    x1 = x.x1;
    x2 = x.x2;
else
    x1 = x(:,1);
    x2 = x(:,2);
end

y = exp(0.7*x1 + 0.3*x2);