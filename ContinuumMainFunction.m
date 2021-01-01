function [E0,ADRatio] = ContinuumMainFunction(country)
%run this file from the command window with one of the following string inputs to reproduce Table 4:
% 'Italy' , 'Germany', 'UK', 'Spain'
global nDays data_cal data ys sigma cf gamma1 tolN
if strcmp(country,'Italy'); X=[-2.00188754602483,1.62595165533504,1.95047276782753,2.29892201682886,-1.61792666323166]; 
elseif strcmp(country,'Germany'); X=[-1.89522721309271,1.67228015713143,1.97912619743166,-0.153459534071443,-2.54704622412673];
elseif strcmp(country,'UK'); X=[-2.13319452757009,1.58964187921160,1.71034216222796,-1.46787516089771,-0.740241543278229];
elseif strcmp(country,'Spain'); [-2.94303913317586,1.78072240689756,1.31624278957644,1.00759547374750,-4.00423791046035];
end

sigma=3; gamma1=0.07; tolN=8; 
[data,ys,cf]=countryData(country); 
data=data/ys/1e8; nDays=numel(data); 
[data_cal, ADRatio]=ddeContinuum(X);
e=data_cal-data(cf:end);
E0=norm(e)/norm(data(cf:end));
E0=100*E0;
end