%function  tau = ID( model, q, qd, qdd, f_ext )

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.

%% Modified version of space motion based Newton Euler with change of variable for efficient calculation purpose
%%
clc,close all,clear all
echo off;
%% Configuration Setup
WriteMatlabFunction=false;%May take time
MassMatrix=false;
IncludeMass=true;
IncludeInertia=true;

ComputeRegressor=true;
InertiaMatrixCalculation=true;
NumericalInitialization=true;
SymbolicInitialization=not(NumericalInitialization);
RunScriptTest=false;
Simplification=false;
UseReplacementVariables=true;
InertiaMatrixandCoriolisCalculation=true;
RunID=num2str(round(abs(randn(1)*1E6)));
%%
Ndof=10;
Nbodies=Ndof;
Njoints=Ndof;
addpath('C:\Users\Soldà Enrico\OneDrive\Desktop\Lauree 2022\I Semester 2023\Thesis\Simscape Example\')
addpath('C:\Users\Soldà Enrico\OneDrive\Desktop\Lauree 2022\I Semester 2023\Thesis\Simscape Example\Tree')

Tree=char(ones(1,Ndof)*'R');
Ref=[1:length(Tree)];%Ref(4)=2;Ref(6)=2;
[DH,q,L]=BuildDHSequence(Tree,Ref);
[DH,q,L]=BuildDHSequenceFrom0aAlphadth(Tree,Ref);
clear q
q=sym('q%d',[Ndof,1],'real');
dq=sym('dq%d',[Ndof,1],'real');
ddq=sym('ddq%d',[Ndof,1],'real');

p=[DH{:,5}]-1;
jointType=[DH{:,7}];
Nbodies=size(DH,1);
Ndof=Nbodies;
Nmasses=1;
SaveAll=false;
coordRange=1:3;
dofRange=1:Ndof;
EnableOffsets=1;
TimeCutOff=500;
%% Definition of auxiliary functions

%%
base_acc = str2sym('[dw0x*1;dw0y*1;dw0z*1;g0x;g0y;g0z]');

syms t
% q=transpose(GetSimVec(Ndof,{'q'},1,'(t)'));
% dq=diff(q,t);
% ddq=diff(dq,t);
%%

T{1}=eye(6);
%%
rotz=@(q)[  cos(q)  sin(q)  0  0  0  0 ;...
    -sin(q)  cos(q)  0  0  0  0 ;...
    0  0  1  0  0  0 ;...
    0  0  0  cos(q)  sin(q)  0 ;...
    0  0  0 -sin(q)  cos(q)  0 ;...
    0  0  0  0  0  1 ];
rotx=@(q)[ 1 0 0 0 0 0;...
    0 cos(q)  sin(q)   0  0  0 ;...
    0 -sin(q)  cos(q)  0  0  0 ;...
    0  0  0  1  0 0;...
    0  0  0  0  cos(q)  sin(q);...
    0  0  0  0 -sin(q)  cos(q);...
    ];
xlt=@(r)[  1     0     0    0  0  0 ;...
    0     1     0    0  0  0 ;...
    0     0     1    0  0  0 ;...
    0     r(3) -r(2) 1  0  0 ;...
    -r(3)  0     r(1) 0  1  0 ;...
    r(2) -r(1)  0    0  0  1 ];

DH2XT=@(Alpha,a,d,theta)[xlt([a; 0; d])*rotx(Alpha)*rotz(theta)];
tilde=@(v)[0 v(3) -v(2); -v(3) 0 v(1);v(2) -v(1) 0];

% rbi = [ I + m, m*tilde(C); m*transpose(tilde(C)), m*eye(3) ];
rbi = @(I,C,m)[ I + m*tilde(C)*transpose(tilde(C)), m*tilde(C); m*transpose(tilde(C)), m*eye(3) ];
Mlink=@(m,I)[I,zeros(3);zeros(3),m*eye(3)];
%%
for k=1:Nbodies
    % eval(['T{' num2str(DH{k,6}) '}=DHProximal(DH{k,1},DH{k,2},DH{k,3},DH{k,4})'])
    T{k}=DH2XT(DH{k,1},DH{k,2},DH{k,3},DH{k,4});
    XL{k}=DH2XT(DH{k,1},-DH{k,2},DH{k,3}*0,DH{k,4}*0);
end
%%
for k=1:Nbodies
    rcm(coordRange,k)=str2sym(['[rcm' num2str(k) 'x' ';rcm' num2str(k) 'y*' num2str(EnableOffsets) ';rcm' num2str(k) 'z*' num2str(EnableOffsets) ';]']);%Center of mass in link frame
    % Lp(coordRange,k)=str2sym(['[L' num2str(k) 'x;L' num2str(k) 'y;L' num2str(k) 'z;]']);%Center of mass in link frame
end
if Nmasses>1
    m=GetSimVec([Ndof,Nmasses],{'m'})*IncludeMass;
else
    m=GetSimVec(Ndof,{'m'})*IncludeMass;

end
for k=1:Nbodies
    Ic{k}=1*((GetSimVecSymmetric([3,3],{['I' num2str(k) '_']})))*IncludeInertia;
    I{k}=rbi(Ic{k},rcm(:,k),m(k));
end;
%% Initialization
v{1}=str2sym('[w0x*0;w0y*0;w0z;v0x;v0y;v0z]');
dv{1}=str2sym('[dw0x;dw0y;dw0z;v0x;v0y;v0z]');
disp('Start forward kinematics')
figure('Name','Processing time forward')
plot(0,0,'sr');hold on;grid on;
ts=tic;
Mp=0;
TAUPRIME=[ ];
KINLIST=[];
WRENCHLIST=[];
%% Forward kinematics

% if InertiaMatrixCalculation
% end
for i = 1:Nbodies
    tb=tic;
    [ XJ, S{i} ] = jcalc( jointType(i), q(i) );
    vJ = S{i}*dq(i);
    if i>1
        Xup{i} = simplify(XJ * XL{i});
    else
        Xup{i} = simplify(XJ *eye(6));
    end
    if p(i) == 0
        v{i} = vJ;
        dv{i} = Xup{i}*(-base_acc) + S{i}*ddq(i);
    else
        if Simplification
            v{i} = simplify(Xup{i}*v{p(i)} + vJ);
        else
            v{i} = Xup{i}*v{p(i)} + vJ;
        end
        dv{i} = Xup{i}*dv{p(i)} + S{i}*ddq(i) + crm(v{i})*vJ;
    end

    % Calculate jacobian and mass matrix components
    % J{i}=jacobian(v{i},dq);
    % Mp=Mp+transpose(J{i})*Mlink(m(i),Ic{i})*J{i};
    %%
    disp(['Body ' num2str(i) ': done.']);
    tproc(i)=toc(tb);
    if or((tproc(i)>TimeCutOff),UseReplacementVariables)
        disp('Change of variable.');
        eval(['syms ' 'w' num2str(i) '_1 w' num2str(i) '_2 w' num2str(i) '_3 v' num2str(i) '_1 v' num2str(i) '_2 v' num2str(i) '_3']);
        eval(['syms ' 'dw' num2str(i) '_1 dw' num2str(i) '_2 dw' num2str(i) '_3 dv' num2str(i) '_1 dv' num2str(i) '_2 dv' num2str(i) '_3']);
        eval(['KINLIST=[KINLIST;[dw' num2str(i) '_1; dw' num2str(i) '_2; dw' num2str(i) '_3]==dv{' num2str(i) '}(1:3);[dv' num2str(i) '_1; dv' num2str(i) '_2; dv' num2str(i) '_3]==dv{' num2str(i) '}(4:6);[w' num2str(i) '_1; w' num2str(i) '_2; w' num2str(i) '_3]==v{' num2str(i) '}(1:3);[v' num2str(i) '_1; v' num2str(i) '_2; v' num2str(i) '_3]==v{' num2str(i) '}(4:6);]' ]);

        % kin_def_expr=['w' num2str(i) '_1==v{' num2str(i) '}(1);w' num2str(i) '_2==v{' num2str(i) '}(2);w' num2str(i) '_3==v{' num2str(i) '}(3);v' num2str(i) '_1==v{' num2str(i) '}(4);v' num2str(i) '_2==v{' num2str(i) '}(5);v' num2str(i) '_3==v{' num2str(i) '}(6);'];
        % eval(kin_def_expr);
        eval(['w' num2str(i) '_1=v{' num2str(i) '}(1);w' num2str(i) '_2=v{' num2str(i) '}(2);w' num2str(i) '_3=v{' num2str(i) '}(3);v' num2str(i) '_1=v{' num2str(i) '}(4);v' num2str(i) '_2=v{' num2str(i) '}(5);v' num2str(i) '_3=v{' num2str(i) '}(6);'])
        eval(['dw' num2str(i) '_1=dv{' num2str(i) '}(1);dw' num2str(i) '_2=dv{' num2str(i) '}(2);dw' num2str(i) '_3=dv{' num2str(i) '}(3);dv' num2str(i) '_1=dv{' num2str(i) '}(4);dv' num2str(i) '_2=dv{' num2str(i) '}(5);dv' num2str(i) '_3=dv{' num2str(i) '}(6);'])
        eval(['w' num2str(i) '=v{' num2str(i) '}(1:3);v' num2str(i) '=v{' num2str(i) '}(4:6);']);
        eval(['dw' num2str(i) '=dv{' num2str(i) '}(1:3);dv' num2str(i) '=dv{' num2str(i) '}(4:6);']);
        % v{i}=str2sym(['[w'  num2str(i) '(1:3);v' num2str(i) ']']);
        % dv{i}=str2sym(['[dw' num2str(i) ';a' num2str(i) ']']);
        %eval(['KINLIST=[KINLIST; w' num2str(i) ';v' num2str(i) '];']);

        v{i}=str2sym(['[w'  num2str(i) '_1;w' num2str(i) '_2;w'  num2str(i) '_3;v' num2str(i) '_1;v' num2str(i) '_2;v' num2str(i) '_3;]']);
        dv{i}=str2sym(['[dw' num2str(i) '_1;dw' num2str(i) '_2;dw' num2str(i) '_3;dv' num2str(i) '_1;dv' num2str(i) '_2;dv' num2str(i) '_3;]']);

    end
    %% Calculation of forcres and torques due to inertia
    f{i} = I{i}*dv{i} + crf(v{i})*I{i}*v{i};
    WrchI{i}= f{i};
    plot(1:i,log(tproc),'-bs');grid on;xlabel('Body #');ylabel('log computation time [s]');
end
tp=toc(ts)
disp('Start backward calculations..')
%% Add external forces
% if nargin == 5
%   f = apply_external_forces( model.parent, Xup, f, f_ext );
% end
% taus=zeros(Ndof,1);
% syms('taus',[Ndof,1]);
taus=[];
plot(0,0,'sg');hold on;grid on;
tp=tic;
for i = Nbodies:-1:1
    tb=tic;
    tau(i,1) = transpose(S{i}) * f{i};%Joint actuation
    eval(['syms taus' num2str(i)]);
    eval(['taus=[taus' num2str(i) '==tau(' num2str(i) ',1); taus]']);
    if p(i) ~= 0
        f{p(i)} = f{p(i)} + transpose(Xup{i})*f{i};

        disp(['Step forces on body ' num2str(i) ': Done.']);
        tprocb(Nbodies-i+1)=toc(tb);
        if or(tprocb(Nbodies-i+1)>TimeCutOff,UseReplacementVariables)
            disp('Change of variable.');
            eval(['syms ' 'm' num2str(i) '_1 m' num2str(i) '_2 m' num2str(i) '_3 f' num2str(i) '_1 f' num2str(i) '_2 f' num2str(i) '_3']);%Declare variables for force-moments vector
            eval(['syms ' 'm' num2str(p(i)) '_1 m' num2str(p(i)) '_2 m' num2str(p(i)) '_3 f' num2str(p(i)) '_1 f' num2str(p(i)) '_2 f' num2str(p(i)) '_3']);%Declare variables for force-moments vector

            eval(['WRENCHLIST=[WRENCHLIST;  m' num2str(p(i)) '_1==f{' num2str(p(i)) '}(1);m' num2str(p(i)) '_2==f{' num2str(p(i)) '}(2);m' num2str(p(i)) '_3==f{' num2str(p(i)) '}(3);f' num2str(p(i)) '_1==f{' num2str(p(i)) '}(4);f' num2str(p(i)) '_2==f{' num2str(p(i)) '}(5);f' num2str(p(i)) '_3==f{' num2str(p(i)) '}(6);]']);

            % eval(['WRENCHLIST=[WRENCHLIST;[m' num2str(i) '_1; m' num2str(i) '_2; m' num2str(i) '_3]==f{' num2str(i) '}(1:3);[f' num2str(i) '_1; f' num2str(i) '_2; f' num2str(i) '_3]==f{' num2str(i) '}(4:6);]']);

            % wrench_def_expr=['m' num2str(p(i)) '_1==f{' num2str(p(i)) '}(1);m' num2str(p(i)) '_2==f{' num2str(p(i)) '}(2);m' num2str(p(i)) '_3==f{' num2str(p(i)) '}(3);f' num2str(p(i)) '_1==f{' num2str(p(i)) '}(4);f' num2str(p(i)) '_2==f{' num2str(p(i)) '}(5);f' num2str(p(i)) '_3==f{' num2str(p(i)) '}(6);'];
            % eval(wrench_def_expr);

            eval(['m' num2str(p(i)) '_1=f{' num2str(p(i)) '}(1);m' num2str(p(i)) '_2=f{' num2str(p(i)) '}(2);m' num2str(p(i)) '_3=f{' num2str(p(i)) '}(3);f' num2str(p(i)) '_1=f{' num2str(p(i)) '}(4);f' num2str(p(i)) '_2=f{' num2str(p(i)) '}(5);f' num2str(p(i)) '_3=f{' num2str(p(i)) '}(6);'])
            % eval(['m' num2str(p(i)) '=f{' num2str(p(i)) '}(1:3);f' num2str(p(i)) '=f{' num2str(p(i)) '}(4:6);']);

            %eval(['WRENCHLIST=[WRENCHLIST; m' num2str(p(i)) ';f' num2str(p(i)) '];']);

            f{p(i)}=str2sym(['[m' num2str(p(i)) '_1;m' num2str(p(i)) '_2;m' num2str(p(i)) '_3;f' num2str(p(i)) '_1;f' num2str(p(i)) '_2;f' num2str(p(i)) '_3];']);
            % f{p(i)}=str2sym(['[m' num2str(p(i)) ';f' num2str(p(i)) '];']);
            %eval(['WRENCHLIST=[WRENCHLIST;[m' num2str(i) '_1; m' num2str(i) '_2; m' num2str(i) '_3]==f{' num2str(i) '}(1:3);[f' num2str(i) '_1; f' num2str(i) '_2; f' num2str(i) '_3]==f{' num2str(i) '}(4:6);]']);

        end
    end
    plot(log(tprocb),'-gs');grid on;xlabel('Body #');ylabel('log computation time [s]');
end
TAUPRIME=[KINLIST;WRENCHLIST;taus];
legend('Forward','Backward');
tbk=toc(tp)
%% Write initialization file section
fileID = fopen(['S' RunID '_' 'Init_' 'Tmp.m'],'w');
if NumericalInitialization
    fprintf(fileID,'%% =========  INTIALIZATION ================ \n');
    fprintf(fileID,'g0x=0;g0y=0;g0z=9.81;   \n');
    fprintf(fileID,'w0x=0;w0y=0;w0z=0;   \n');
    fprintf(fileID,'dw0x=0;dw0y=0;dw0z=0;   \n');
    fprintf(fileID,'\n');
elseif SymbolicInitialization
    fprintf(fileID,'syms g0x g0y g0z  \n');
    fprintf(fileID,'syms w0x w0y w0z   \n');
    fprintf(fileID,'syms dw0x dw0y dw0z   \n');
end

for k=1:Ndof
    if NumericalInitialization
        fprintf(fileID,'m%d=1; Alpha%d=0; a%d=1; d%d=1; q%d=rand(1); \n',k,k,k,k,k);
        fprintf(fileID,'dq%d=rand(1); ddq%d=rand(1); \n',k,k);
        fprintf(fileID,'I%d_11=1; I%d_22=1; I%d_33=1; I%d_12=1; I%d_13=1; I%d_23=1;  \n',k,k,k,k,k,k);
        fprintf(fileID,'rcm%dx=rand(1); rcm%dy=rand(1); rcm%dz=rand(1);  \n',k,k,k);
        fprintf(fileID,'\n');
    elseif SymbolicInitialization
        fprintf(fileID,'syms m%d Alpha%d a%d d%d q%d \n',k,k,k,k,k);
        fprintf(fileID,'syms dq%d ddq%d \n',k,k);
        fprintf(fileID,'syms I%d_11 I%d_22 I%d_33 I%d_12 I%d_13 I%d_23  \n',k,k,k,k,k,k);
        fprintf(fileID,'syms rcm%dx rcm%dy rcm%dz  \n',k,k,k);
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
fprintf('\n Initialization part written. \n')
%% write equations on file
fileID = fopen(['S' RunID '_' 'Tmp.m'],'a');
fprintf(fileID,'%% =========  EXECUTION SECTION =========== \n');
% Convert to writable and executable string
KINLIST4print=strrep(string(KINLIST),'==','=');
WRENCHLIST4print=strrep(string(WRENCHLIST),'==','=');
tau4print= strrep(string(taus),'==','=');

fprintf(fileID,'\n %s; \n %s \n',KINLIST4print,'%Comment: End kinematics \n');

% for k=1:length(WRENCHLIST)
%     fprintf(fileID,'%s \n',strrep(string(WRENCHLIST(k)),'==','='));
% end
fprintf(fileID,'\n %s; \n %s \n',WRENCHLIST4print,'%Comment: End whrench  \n');
fprintf(fileID,'\n %s; \n  %s \n',tau4print,'%Comment: End tau \n');
fclose(fileID);
fprintf('\n Execution part written. \n')

%% Check script for running
if and(or(NumericalInitialization,Ndof<8),RunScriptTest)%avoid execution in symbolic if long expressions
fprintf('Excution tryal... \n')
run(['S' RunID '_Init_' 'Tmp.m']);
run(['S' RunID '_' 'Tmp.m']);
fprintf('Done. \n');
pause(5);
else
    warning('Skipped evaluation of script')
end
%% Post processing
% disp('Build model structure')
% model.TAUPRIME=TAUPRIME;
% model.tau=tau;
% model.vars=symvar(TAUPRIME);
% fun=matlabFunction(TAUPRIME,'File','Dummy');
% disp('Done.');

%% Save and report
if SaveAll
    disp('Saving..');
    SessionName=['RNEA_Space_' RunID '_' num2str(Ndof)];
    save(SessionName);
    disp('Done.');
end

%% Inertia matrix calculation
if InertiaMatrixandCoriolisCalculation
    disp('Calculating algebrically mass matrix and coriolis terms...')
    M=jacobian(tau,ddq);
    C=jacobian(tau-M*ddq,dq);
end

if WriteMatlabFunction
    EOMhandle=matlabFunction(tau,'File',[SessionName '_EOM' ]);
end
%% function

function [Xj,S]= jcalc(jtype,q)
switch jtype
    case 0	% revolute Z axis
        Xj = rotz(q);
        S = [0;0;1;0;0;0];
    case 1			% prismatic Z axis
        Xj = xlt([0 0 q]);
        S = [0;0;0;0;0;1];
    case 2				% helical (Z axis)
        Xj = rotz(q) * xlt([0 0 q*jtyp.pars.pitch]);
        S = [0;0;1;0;0;jtyp.pars.pitch];
    otherwise
        error( 'unrecognised joint code ''%s''', code );
end
end