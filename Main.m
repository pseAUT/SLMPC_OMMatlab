
close all
clear all
clc
               


global Ac Bc Cc Dc X_S derx_S Y_S u_lin
global NH NC Ts Xhat nus nuo Yhat hsp kc Adis Bdis Kdis Jdis QQ1 QQ2 u_pre

hsp=[0.56;0.52;0.36];

QQ1=diag([1;1;1]);

QQ2=diag([18;18;18]);

u_pre=[0;0;0];
import OMMatlab.*



omc = OMMatlab;
omc.ModelicaSystem("tank.mo","tank");

omc1 = OMMatlab;
omc1.ModelicaSystem("tank.mo","tank");

%%
% Initial flow Inputs
I1=0;    
I2=0;    
I3=0;

u_lin=[I1;I2;I3];

%% initial states

IC=[0;0;0];
    
%% initial parameter names for initialization
State_names=["h1_In";"h2_In";"h3_In"];

%% state fields, to summon from result files
State_Fields=  ["h1";"h2";"h3"];

%% output variables (coul be summoned by get methods too)
Output_Fields=["qo1";"qo2";"qo3"];
%% 

nus= length(State_Fields);  %number of states

nuo= length(Output_Fields); %number of outputs
%%
%%%%%%%%%%%%%%%%%%%%%%%%%Part by Part%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simulation setup data
Ts=0.1;   % sample time
Tsim=50;  % simulation span
x0=IC; 
n=length(IC);

%% Primary initialization
    for j=1:n
   omc.setParameters(State_names(j)+"="+mat2str(x0(j),15))
    end

%% primary simulation to get all variable values at t=0
omc.setSimulationOptions(["startTime="+num2str(0),"stopTime="+num2str(0),"tolerance=1e-6"]);
omc.simulate();

% pre-specifying Observed outputs
Yhat= zeros(nuo,1);
Y= zeros(nuo,1);

% get initial outputs from primary initialzation
for j1=1:length(Output_Fields)

        tmp=cell2mat(omc.getSolutions(Output_Fields(j1)));
        Y(j1)=tmp(1);    %%%% set initial observed outputs

end
Yhat=Y;

Xhat=[0.1;0.1;0.1];  % initial estimates
Xhat= reshape(Xhat,length(Xhat),1);


%% pre-locating discrete input signal
U=I1*ones(3,ceil(Tsim/Ts));
%%


%% dynamic matrices to store observed data
X_ob=[Xhat];
Y_ob=[Yhat];
T_ob=[0];
X_Real=[x0];
Y_Real=[Y];
T=[0];


%% Noise covariances

% Q=5e-4*eye(nus);
% R=1e-3*eye(nuo);
 Q=1e-4*eye(nus);
 R=2e-4*eye(nuo);

FS_Q=Q(1)*ones(nus,1);



P0=0.1*eye(nus);  % initial State covariance matrix

P=P0;

% Discrete derivative holder

Ad_Pr=zeros(nus,nus);



    Y_S=zeros(nuo,1);


    X_S=zeros(nus,1);

  


    derx_S=zeros(nus,1);

   
%% Set MPC Settings

NH=15;
NC=10;

%%

LB=0*ones(3*NC,1);
UB=0.3*ones(3*NC,1);
U_guess=(LB+UB)/2;
opt_options = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',1e-6, ...
    'MaxIterations',500000,'MaxFunctionEvaluations',10000000,'OptimalityTolerance',1e-6,'StepTolerance',1e-6,...
    'FunctionTolerance',1e-6);

Success_flag=ones(1,ceil(Tsim/Ts));

Dicr_Flag=zeros(1,ceil(Tsim/Ts));
%% do simulation+observation+control:


for i=1:ceil(Tsim/Ts)
tic



    %% initialize the plant model
for j=1:n
omc.setParameters(State_names(j)+"="+mat2str(x0(j),15))
end




    for j=1:n
omc1.setParameters(State_names(j)+"="+mat2str(Xhat(j),15))
    end
%% setting input for omc 
omc.setInputs("qi1="+mat2str(u_lin(1),12));
omc.setInputs("qi2="+mat2str(u_lin(2),12));
omc.setInputs("qi3="+mat2str(u_lin(3),12));

%% setting input for omc1 

omc1.setInputs("qi1="+mat2str(u_lin(1),12));
omc1.setInputs("qi2="+mat2str(u_lin(2),12));
omc1.setInputs("qi3="+mat2str(u_lin(3),12));



%% Linearize 
    omc1.setLinearizationOptions(["startTime="+num2str(0),"stopTime="+num2str(0)]);
    Coeff= omc1.linearize();
    %% Matrices

    Ac=(cell2mat(Coeff(1)));
    Bc=(cell2mat(Coeff(2)));

    Cc=(cell2mat(Coeff(3)));
    Dc=(cell2mat(Coeff(4)));



    %% state and output names 
    Linear_states=omc1.getLinearStates();
    Linear_outputs=omc1.getLinearOutputs();

    %% get Linearization coefficients
    omc1.setSimulationOptions(["startTime=0","stopTime="+num2str(0)]);
    omc1.simulate();

    Y_S=zeros(nuo,1);

    for l=1:1:nuo
    temp=omc1.getSolutions(Output_Fields(l));
    tpY=cell2mat(temp);
    Y_S(l)=tpY(end);
    end
    


    X_S=zeros(nus,1);

    for l=1:1:nus
    temp=omc1.getSolutions(State_Fields(l));
    tpS=cell2mat(temp);
    X_S(l)=tpS(end);
    end



    derx_S=zeros(nus,1);

    for l=1:1:nus
    temp="der("+string(State_Fields(l))+")";
    %der_Fields(l)=temp;
    temp1=omc1.getSolutions(temp);
    derx=cell2mat(temp1);
    derx_S(l)=derx(end);
    end

     kc=derx_S-Ac*X_S-Bc*u_lin;
    Jdis=Y_S-Cc*X_S-Dc*u_lin;

    [Adis,Bdis,Kdis] = discretize(Ac,Bc,kc,Ts);

    %% MPC Optimizer
   [p,fval,exitflag,output] = fmincon (@OBJt,U_guess,[],[],[],[],LB,UB,@nlcondi,opt_options);

    if exitflag<=0
        warning("The optimizer failed at ST= "+num2str(i))
        U(:,i)=U(:,i-1);
        Success_flag(i)=0;
    end
   % U_guess=p;

    U(1,i)=p(1);
    U(2,i)=p(NC+1);
    U(3,i)=p(2*NC+1);
    u_pre=U(:,i);

    %%IIC=x0(1:end-1);
    u_lin=U(:,i);
    u=U(:,i);
    %% setting input for omc 
omc.setInputs("qi1="+mat2str(u(1),12));
omc.setInputs("qi2="+mat2str(u(2),12));
omc.setInputs("qi3="+mat2str(u(3),12));

%% setting input for omc1 

omc1.setInputs("qi1="+mat2str(u(1),12));
omc1.setInputs("qi2="+mat2str(u(2),12));
omc1.setInputs("qi3="+mat2str(u(3),12));



    omc1.setLinearizationOptions(["startTime="+num2str(0),"stopTime="+num2str(0)]);
    Coeff= omc1.linearize();
    %% Matrices

    Ac=(cell2mat(Coeff(1)));
    %% Discrete derivative
    Ad=expm(Ac*Ts);



    %% move plant model to the new ST
    omc.setSimulationOptions(["startTime="+num2str((i-1)*Ts),"stopTime="+num2str((i)*Ts),"tolerance=1e-6"]);
    omc.simulate();

%%%%%get States  @ new ST

    
    for j=1:nus
        
    tmp1=cell2mat(omc.getSolutions(State_Fields(j)));
    tmp2(j)=tmp1(end);
    end
    

    x0=tmp2'+[sqrt(FS_Q).*randn(nus,1)];

    Puritydp_en=tmp2(end);  % Puritydp at the end of sample time


    tmp=cell2mat(omc.getSolutions("time"));
    T=[T,tmp(end)];

    


    %% Measurement
for j1=1:length(Output_Fields)

        tmp=cell2mat(omc.getSolutions(Linear_outputs(j1)));
        Y(j1)=tmp(end);    %%%% initial outputs
end

Y=Y+ sqrt(R(1))* randn(nuo,1);  % adding noise to measurment

 %% Prediction

 xh0=Xhat;
 

   for j=1:n
    omc1.setParameters(State_names(j)+"="+mat2str(xh0(j),15))
   end
   omc1.setSimulationOptions(["startTime="+num2str((i-1)*Ts),"stopTime="+num2str((i)*Ts),"tolerance=1e-6"]);

    omc1.simulate();

    tmp4=zeros(nus,1);

    for j=1:nus
        
    tmp3=cell2mat(omc1.getSolutions(State_Fields(j)));
    tmp4(j)=tmp3(end);
    end
   
    X_pred= tmp4;

    %X_pred=reshape(X_pred,nus,1);%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_pred=Ad*Xhat + Bd* U(:,i) +Ed;
 P_pred=Ad*P*Ad'+ Q;


% xh0=X_pred;

 
    for j=1:n
    omc1.setParameters(State_names(j)+"="+mat2str(X_pred(j),15))
   end

   omc1.setSimulationOptions(["startTime="+num2str((i)*Ts),"stopTime="+num2str((i)*Ts),"tolerance=1e-6"]);

    omc1.simulate();

    Y_pred=zeros(nuo,1);


    for j1=1:nuo

        tmp=cell2mat(omc1.getSolutions(Output_Fields(j1)));
        Y_pred(j1)=tmp(end);   
    end


      omc1.setLinearizationOptions(["startTime="+num2str((i)*Ts),"stopTime="+num2str((i)*Ts),"tolerance=1e-6"]);

     Coeff= omc1.linearize();
    
     Cc1=cell2mat(Coeff(3));



 %Y_pred=Cc*X_pred + Dc*U(:,i)+Fc;
%% Correction
S=Cc1*P_pred*Cc1'+R;
%K=(P_pred*Cc')*inv(S);  % worked
K=(P_pred*Cc1')/S; % pinv is to handle Badly Scaled S effect
Xhat= X_pred+K*(Y-Y_pred);
P=P_pred-K*S*K';



    for j=1:n
    omc1.setParameters(State_names(j)+"="+mat2str(Xhat(j),15))
    end

   omc1.setSimulationOptions(["startTime="+num2str((i)*Ts),"stopTime="+num2str((i)*Ts),"tolerance=1e-6"]);

    omc1.simulate();

    Yhat=zeros(nuo,1);

    for j1=1:nuo
        tmp=cell2mat(omc1.getSolutions(Output_Fields(j1)));
        Yhat(j1)=tmp(end);   
    end


X_ob=[X_ob,Xhat];
Y_ob=[Y_ob,Yhat];
T_ob=[T_ob,i*Ts];
X_Real=[X_Real,x0];
Y_Real=[Y_Real,Y];

disp("Ending ST=  "+num2str(i))

% %%
save("ws_"+num2str(i),'T','T_ob','Y_ob','X_ob','U','Y_Real','X_Real','Dicr_Flag','x0','P','Ts','Tsim')
toc
end
figure(1)

subplot (3,1,1)
stairs(0:Ts:Tsim,[U(1,:) U(1,end)],'LineWidth',6)
set(gca,'fontsize',18)
title('q_{in1}','FontSize',25)
xlabel('Time (min)')
ylabel('  Flow (m^3/min)')
grid on

subplot (3,1,2)
stairs(0:Ts:Tsim,[U(2,:) U(2,end)],'LineWidth',6)
set(gca,'fontsize',18)
title('q_{in2}','FontSize',25)
xlabel('Time (min)')
ylabel('  Flow (m^3/min)')
grid on


subplot (3,1,3)
stairs(0:Ts:Tsim,[U(3,:) U(3,end)],'LineWidth',6)
set(gca,'fontsize',18)
title('q_{in3}','FontSize',25)
xlabel('Time (min)')
ylabel('  Flow (m^3/min)')
grid on


figure(2)
subplot (3,1,1)
plot(T,X_Real(1,:),T_ob,X_ob(1,:),'--','LineWidth',6)
set(gca,'fontsize',18)
legend('Real-Value','Estimated-Value')
title('h_1','FontSize',25)
xlabel('Time (min)')
ylabel('Level(m)')
grid on


subplot (3,1,2)
plot(T,X_Real(2,:),T_ob,X_ob(2,:),'--','LineWidth',6)
set(gca,'fontsize',18)
legend('Real-Value','Estimated-Value')
title('h_2','FontSize',25)

xlabel('Time (min)')
ylabel('Level(m)')

grid on



subplot (3,1,3)
plot(T,X_Real(3,:),T_ob,X_ob(3,:),'--','LineWidth',6)
set(gca,'fontsize',18)
legend('Real-Value','Estimated-Value')
title('h_3','FontSize',25)
xlabel('Time (min)')
ylabel('Level(m)')
grid on


figure(3)
subplot (3,1,1)
plot(T,Y_Real(1,:),'LineWidth',6)
set(gca,'fontsize',18)
legend('Measured-Value')
title('q_{o1}','FontSize',25)
xlabel('Time (min)')
ylabel('  Flow (m^3/min)')
grid on 


subplot (3,1,2)
plot(T,Y_Real(2,:),'LineWidth',6)
set(gca,'fontsize',18)
legend('Measured-Value')
grid on
title('q_{o2}','FontSize',25)
xlabel('Time (min)')
ylabel('  Flow (m^3/min)')

subplot (3,1,3)
plot(T,Y_Real(3,:),'LineWidth',6)
set(gca,'fontsize',18)
legend('Measured-Value')
title('q_{o3}','FontSize',25)
xlabel('Time (min)')
ylabel('  Flow (m^3/min)')
grid on
