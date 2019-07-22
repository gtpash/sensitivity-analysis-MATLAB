%
%  Code for the Morris algorithm to construct global sensitivity indices.
%

  clear all

%
% Use a stepsize of ell = 1/100 and r evaluations.
%

  r = 200;
  ell = 100;

  val1 = 1;
  val2 = 1;
  
  rng(5);

  for i = 1:r
    Num = rand(1,2);
    for j = 1:ell-2
      v1 = abs(j/ell - Num(1));
      v2 = abs(j/ell - Num(2));
      if v1 < val1
         q1 = j/ell;
         val1 = v1;
      end
      if v2 < val2
         q2 = j/ell;
         val2 = v2;
      end
    end
    q = [q1 q2];
    Q(i,:) = q;
    Delta = ell/(2*(ell-1));
    
    num = rand(1);
    if 0 < num & num < 0.5
      P = [1 0;0 1];
    else
      P = [0 1; 1 0];
    end
%     P
    
    for j=1:2
      num = rand(1);
      if num < 0.5
        D(j,j) = -1;
      else
        D(j,j) = 1;
      end
    end
%     D
    J1 = ones(3,1);
    J3 = ones(3,2);
    B = [0 0;1 0;1 1];
    Bs = (J1*q + (Delta/2)*((2*B - J3)*D + J3))*P;
  
    for j=1:3
      params = Bs(j,:);
      k = params(1);
      m = params(2);
      K = sqrt(k/m);
%      y(j) = cos(K*pi/2);
      y(j) = 1/K;
    end
 
    for j=1:2
      if Bs(j+1,1)~=Bs(j,1)
        d1(i) = (y(j+1) - y(j))/Delta;
      elseif Bs(j+1,2)~=Bs(j,2)
        d2(i) = (y(j+1) - y(j))/Delta;
      end
    end
  end

%
% Construct the Morris indices.
%

  mu1 = (1/r)*sum(d1);
  mu2 = (1/r)*sum(d2);

  mu1s = (1/r)*sum(abs(d1));
  mu2s = (1/r)*sum(abs(d2));
  mus = [mu1s mu2s]

  sigma1 = (1/(r-1))*sum((d1-mu1).^2);
  sigma2 = (1/(r-1))*sum((d2-mu2).^2);
  sigma = [sigma1 sigma2]