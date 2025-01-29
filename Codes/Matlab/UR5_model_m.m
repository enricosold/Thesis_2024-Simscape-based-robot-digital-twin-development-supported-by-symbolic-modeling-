%%
Ndof=9;
Run_ID=num2str(round(randn(1)*1E7));
%%
EnableFriction=true;
g=-9.81;

%%Geometrical parameters
lf1=30;
df=6;
lf2=lf1;
lf3=lf1;

L=0.5;
d=0.05;
h0=0.05;
r=0.05;
h=0.1;
wl=h;
wl2=0.15;
w1=0;

L1=0.16;
L2=.45;
L3=0.4;
L4=0.15;
L5=0.15;
L6=0.15;
L7=0.15;

%% Inertial parameters
m0=1e-15;%Joint mass
m=1;
Cm=[0 0 0 0 0 0 0 0;...
   0 0 0 0 0 0 0 0;...
   0 0 0 0 0 0 0 0];
%% UR5 link masses and geometry
% L=[];
% d1=0.1625;
% a2=-0.425;
% a3=-0.3922;
% d4=	0.1333;
% d5=0.0997;
% d6=0.0996;

rcm1=[0, -0.02561, 0.00193];
rcm2=[0.2125, 0, 0.11336];
rcm3=[0.15, 0.0, 0.0265];
rcm4=[0, -0.0018, 0.01634];
rcm5=[0, 0.0018,0.01634];	
rcm6=[0, 0, -0.001159];

ma1=3.7%Large actuators at the beginning of chain
ma2=1.2;%Small actuators at the end of chain
ma2=1.0;
ma2=1.25;
Mtot=([ 3.761 	8.058 	2.846 1.37 	1.3 0.365  0 0 0 0 0 0]+0)*m;%Total link masses
M=[0 1.28 1.28 0 0 0.2];%Link masses (neglected firt, fourth fifth and sixth)
% M=[0 8 2.8 0 0 0];%Link masses (neglected firt, fourth fifth and sixth)

Mact=([ ma1 ma1 ma1 ma2 ma2 ma2])*1;%Actuator mass
mp=0.05;
%% UR5 model has negligible inertia
I1z=1*0;
I2z=0.6*0;
I3z=0.5*0;
I4z=0.3*0;
I5z=0.114*1;
I6z=0.100*1;

I1=[[0 0   I1z] [0 0 0]];
I2=[[0 0   I2z] [0 0 0]];
I3=[[0 0   I3z] [0 0 0]];
I4=[[0 0   I4z] [0 0 0]];
I5=[[0 0   I5z] [0 0 0]];
I6=[[0 0   I6z] [0 0 0]];
% I7=[[0 0 0] [0 0 0]];

%% Friction parameters
EnableFrictionInitialValues=false;
Fc1=0.1*EnableFrictionInitialValues;
Fc2=0.1*EnableFrictionInitialValues;
Fc3=0.1*EnableFrictionInitialValues;
Fc4=0.1*EnableFrictionInitialValues;
Fc5=0.1*EnableFrictionInitialValues;
Fc6=0.1*EnableFrictionInitialValues;

Fv1=0.1*EnableFrictionInitialValues;
Fv2=0.1*EnableFrictionInitialValues;
Fv3=0.1*EnableFrictionInitialValues;
Fv4=0.1*EnableFrictionInitialValues;
Fv5=0.1*EnableFrictionInitialValues;
Fv6=0.1*EnableFrictionInitialValues;

Fvs1=0.0*EnableFrictionInitialValues;
Fvs2=0.0*EnableFrictionInitialValues;
Fvs3=0.0*EnableFrictionInitialValues;
Fvs4=0.0*EnableFrictionInitialValues;
Fvs5=0.0*EnableFrictionInitialValues;
Fvs6=0.0*EnableFrictionInitialValues;


Fvc1=0.0*EnableFrictionInitialValues;
Fvc2=0.0*EnableFrictionInitialValues;
Fvc3=0.0*EnableFrictionInitialValues;
Fvc4=0.0*EnableFrictionInitialValues;
Fvc5=0.0*EnableFrictionInitialValues;
Fvc6=0.0*EnableFrictionInitialValues;


%% Jont laws
Kc=1000*0;
Dc=10*0;
run('UR5_setup')
%% Solid colors
LinkColor =[0 0.8 0.8];
Jointcolor=[1 0 0.2];
Color=[0.0 0.5 0.5];
Opacity=1;
%% Initial conditions
q0=[0;0;0;0;0;0];dq0=[0;0;0;0;0;0];ddq0=[0;0;0;0;0;0];
%% Build input for simulation
syms C1 C2 w0 t
qe = C1.*sin(w0*t)+C2.*sin(w0*t);
dqe=diff(qe,t);
ddqe=diff(dqe,t);
% Numerical istantiation
C1=([1:Ndof]);
C2=randn([1,Ndof]);
w0=[1:Ndof]*0.1;
t=0:0.1:10;t=transpose(t);
qh=(eval(qe));
dqh=(eval(dqe));
ddqh=(eval(ddqe));
figure('Name','Inputs')
plot(t,qh,t,dqh,t,ddqh);grid on;xlabel('t [s]')
%%
open_system('GeneralRobot.slx')
StopTime=max(t);
RobotName='UR5_model_simulink';
% RobotName='UR5_model_simulink_force_driven';%Force driven model
mdl=RobotName;
open_system(RobotName);

%%
simIn = Simulink.SimulationInput(mdl);
simIn = setModelParameter(simIn,"Solver","daessc",...
    "StopTime","100");
out=sim(simIn);%Simulate model using assigned parameters
%% Save data and filter
save([RobotName '_store_' RunID])
%% Run filter in postprocessing
Np=21;
P=diag(1E8*ones(1,Np));
for k=1:length(t)
    obs=reshape(out.tau(3*[1:6],1,k),[Ndof,1]);
    
end