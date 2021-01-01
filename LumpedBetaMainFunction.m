function [E0,ADRatio,R0] = LumpedBetaMainFunction(country) 
%run this file from the command window with one of the following string inputs to reproduce Table 3:
% 'Italy' , 'Germany', 'UK', 'Spain'
global nDays data_cal data ys sigma cf gamma1 tolN
if strcmp(country,'Italy'); X=[-1.70077176705963,2.11973581888038,0.395719351626882,-0.721260697939164]; 
elseif strcmp(country,'Germany'); X=[-1.56197393248550,2.23182335036030,-0.119353967471291,-0.651883083842769];
elseif strcmp(country,'UK'); X=[-1.81013042141519,2.07079833806585,0.385162980548601,-0.739544469042364];
elseif strcmp(country,'Spain'); X= [-1.72308172558173,1.98492630495309,-0.343679650958787,-0.302993722634888];
end

sigma=3; gamma1=0.07; tolN=8; 
beta=exp(X(1)); p=exp(-X(2)^2); tau=1+13*exp(-X(3)^2); pbar=p*(1-exp(-gamma1*tau));
[data,ys,cf]=countryData(country); 
data=data/ys/1e8; nDays=numel(data); 
[data_cal,ADRatio]=ddeLumped(X);
e=data_cal-data(cf:end);
E0=norm(e)/norm(data(cf:end));
E0=100*E0; R0=beta*(1-pbar)/gamma1;
end