%one-step discrete integration
function [x,y]=predModel(x0,u)
global Adis Bdis Cc Dc Kdis Jdis

x=Adis*x0 + Bdis*u+Kdis;

y=Cc*x+Dc*u+Jdis;
end
