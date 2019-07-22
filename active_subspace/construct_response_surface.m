function f_approx_eval = construct_response_surface(qi, y, num)
%% Construct response surface
% Using regression analysis construct a response surface r(y), using q_i,
% such that q_i ~ r(y_i).
poly_order = 3;
f_approx = polyfitn(y',qi,poly_order);
y_eval = linspace(min(y),max(y),num);
f_approx_eval = polyvaln(f_approx,y_eval);
end