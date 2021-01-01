function [Z1,Z2]=ddeContinuum(X)
global tolN data nDays sigma tau gamma1 cf 
a=exp(X(1)); m=2.1+8*exp(-X(2)^2); p=exp(-X(3)^2); tau=1+13*exp(-X(4)^2); f0=1e-2*exp(X(5)); 
nu=sigma+tau; pbar=p*exp(-gamma1*tau);
op=ddeset('reltol',10^(-tolN),'abstol',10^(-tolN));
lags=[sigma nu]; 

tspan=[cf nDays];
sol=dde23(@ddeSet,lags,@hist,tspan,op);
h=sol.y(2,:); t=sol.x;
Z1=interp1(t,h,cf:nDays);

tspan=[cf 2*nDays];
sol=dde23(@ddeSet,lags,@hist,tspan,op);
f=sol.y(1,:); h=sol.y(2,:); S=(a*f).^(m-1)*(m-1).*exp(a*f).*gamma_incomplete(a*f,1-m);
Z2=(1-S(end))/h(end);

function z = ddeSet(~,y,Z)
ylag1=Z(1,1); ylag2=Z(1,2);
z=zeros(2,1); 
    function zG=Gfun(yArg)
             g2A=gamma_incomplete(abs(a*yArg),2-m);
             zG=yArg^(m-1)*a^m*exp(a*yArg)*g2A+yArg^(-2+m)*a^(m-1)*exp(a*yArg)*g2A*m-a^(m-1)*yArg^(-2+m)*exp(a*yArg)*g2A-a;
    end

z(1)=-Gfun(ylag1)+pbar*Gfun(ylag2)-gamma1*y(1)+(1-pbar)*a/(m-2);
z(2)=pbar*z(1)*Gfun(y(1));
end

  function z = hist(~)
      z=[f0; data(cf)];
  end
end