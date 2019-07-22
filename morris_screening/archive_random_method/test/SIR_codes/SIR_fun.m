%
%            function SIR_fun
%

  function  Y = SIR_fun(time,params)

  global N Y0

  ode_options = odeset('RelTol',1e-2);
  [t,Y] = ode45(@SIR_rhs,time,Y0,ode_options,params);
