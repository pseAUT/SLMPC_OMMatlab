function [c,ceq] = nlcondi(uinp)
global Ac Bc Cc Dc X_S derx_S Y_S u_lin
global NH NC Ts Xhat Yhat nus nuo hsp md Q1 Q2


p=[reshape(uinp(1:NC),1,NC);reshape(uinp(NC+1:2*NC),1,NC);reshape(uinp(2*NC+1:end),1,NC)];

p=[p,p(:,end).*ones(3,NH-NC)];
x0=Xhat;

epsi=1e-3;

h_limit=1;


tc(1)=0;
tc(2)=0;
% tc(3)=0;
% tc(4)=0;
% tc(5)=0;


for i=1:NH
[x,y]=predModel(x0,p(:,i));

x0=x;


tc(1)=tc(1)+max(0,-y(1));
tc(2)=tc(2)+max(0,-y(2));

% tc(3)=tc(3)+max(0,x(1)-h_limit);
% tc(4)=tc(4)+max(0,x(2)-h_limit);
% tc(5)=tc(5)+max(0,x(3)-h_limit);

end
c(1)=tc(1)-epsi;
c(2)=tc(2)-epsi;
% 
% c(3)=tc(3)-epsi;
% c(4)=tc(4)-epsi;
% c(5)=tc(5)-epsi;

ceq=[];

end