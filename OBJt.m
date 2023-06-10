function J = OBJt(uinp)

global Ac Bc Cc Dc X_S derx_S Y_S u_lin
global NH NC Ts Xhat Yhat nus nuo hsp md QQ1 QQ2 u_pre


p=[reshape(uinp(1:NC),1,NC);reshape(uinp(NC+1:2*NC),1,NC);reshape(uinp(2*NC+1:end),1,NC)];

p=[p,p(:,end).*ones(3,NH-NC)];


Memorizedp_added=[u_pre,p];

T=[];
X=[];
Y=[];

 


size_in=size(p);

x0=Xhat;
%y0=(Cc*(x0-X_S)+Dc*(p(:,1)-u_lin)+Y_S);

J=0;

for i=1:NH
[x,~]=predModel(x0,p(:,i));

x0=x;
%%du= p(:,i+1)-p(:,i);
%J= J+(x(1)-hsp(1))^2+(x(2)-hsp(2))^2+(x(3)-hsp(3))^2;
du= p(:,i)-Memorizedp_added(:,i);
J= J+(x-hsp)'*QQ1*(x-hsp)+du'*QQ2*du;


end

end