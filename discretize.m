
function [A,B,K] = discretize(Ac,Bc,Kc,Ts)
% For theoritical proof, see Zhakatayev et.al 2017 
    
    A=expm(Ac*Ts);    
    sz=size(A);
    K=Ac\(A-eye(sz))*Kc;              
    B=Ac\(A-eye(sz))*Bc;
end