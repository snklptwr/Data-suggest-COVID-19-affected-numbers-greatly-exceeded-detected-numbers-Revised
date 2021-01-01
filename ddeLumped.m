function [Z1,Z2]=ddeLumped(X)
global nDays data sigma gamma1 tolN cf
beta=exp(X(1)); p=exp(-X(2)^2); tau=1+13*exp(-X(3)^2); in=exp(X(4));
nu=sigma+tau; pbar=p*exp(-gamma1*tau); 
lags=[sigma nu]; op=ddeset('reltol',10^(-tolN),'abstol',10^(-tolN));

tspan=[cf nDays];
sol=dde23(@ddeSet,lags,@hist,tspan,op);
hm=sol.y(2,:); t=sol.x;
Z1=interp1(t,hm,cf:nDays);

tspan=[cf 2*nDays];
sol=dde23(@ddeSet,lags,@hist,tspan,op);
V=sol.y(1,:); hm=sol.y(2,:); S=exp(-beta*V);
Z2=(1-S(end))/hm(end);

function z = ddeSet(~,y,Z)
ylag1=Z(1,1); ylag2=Z(1,2); 
z=zeros(2,1);
z(1)=pbar*exp(-beta*ylag2)-exp(-beta*ylag1)-gamma1*y(1)+1-pbar;
z(2)=beta*pbar*z(1)*exp(-beta*y(1));
end

    function z = hist(~)
        z=[in; data(cf)];
    end
end