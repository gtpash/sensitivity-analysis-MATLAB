
%
%          SIR_rhs
%
  function dy = SIR_rhs(t,y,params);

  theta = params.theta;
  N = params.N;

  gamma = theta.gamma;
  r = theta.r;
  delta = theta.delta;
  k = theta.k;

%
% Exploit the fact that R(t) = N - S(t) - I(t).
%

  dy = [delta*N - delta*y(1) - gamma*k*y(2)*y(1);
        gamma*k*y(2)*y(1) - (r + delta)*y(2)];
