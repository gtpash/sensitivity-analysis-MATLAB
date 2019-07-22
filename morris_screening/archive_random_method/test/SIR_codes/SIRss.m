%
%            function SIRss
%

  function  ss = SIRss(params,data)

  global N Y0

  t_data = data.xdata;
  ydata = data.ydata;

  ode_options = odeset('RelTol',1e-2);
  [t,Y] = ode45(@SIR_rhs,t_data,Y0,ode_options,params);

  ss = sum((Y(:,2) - ydata).^2);