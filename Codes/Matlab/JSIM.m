%% JSIM algorithm for mass matrix calculation
clc,clear all,close all;
profile on

RunID=[num2str(round(randn(1)*1E8))];
diary([RunID '_diary_log'])
diary on
Ndof=4;
DHtype=1;%1: Full DH matrix, 2: simplified DH matrix, 3:Custom
EnableOffsets=0;
Tree=char(ones(1,Ndof)*'R');

Ref=[1:length(Tree)];%Ref(4)=2;Ref(6)=2;%Tree structure
[DH,q,L]=BuildDHSequenceFrom0aAlphadth(Tree,Ref);%less than 9
%[DH,q,L]=BuildDHSequence(Tree,Ref);
for k=1:Ndof,DH{k,3}=0; end;
Nbodies=Ndof;
Nmasses=3;
SimplificationTime=1;
SaveResults=true;
coordRange=1:3;
vecRange=1:6;%Used in space motion
dofRange=1:Ndof;
q=GetSimVec(Ndof,{'q'});
m=GetSimVec(Ndof,{'m'});

%% Kinematics definition
tilde=@(v)[0 v(3) -v(2); -v(3) 0 v(1);v(2) -v(1) 0];

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

rotx =@(a) [ 1  0  0  0  0  0 ;
      0  cos(a)  sin(a)  0  0  0 ;
      0 -sin(a)  cos(a)  0  0  0 ;
      0  0  0  1  0  0 ;
      0  0  0  0  cos(a)  sin(a) ;
      0  0  0  0 -sin(a)  cos(a)];

DH2XT=@(Alpha,a,d,theta)[xlt([a; 0; d])*rotx(Alpha)*rotz(theta)];
  % rbi = [ I + m, m*tilde(C); m*transpose(tilde(C)), m*eye(3) ];
  rbi = @(I,C,m)[ I + m*tilde(C)*transpose(tilde(C)), m*tilde(C); m*transpose(tilde(C)), m*eye(3) ];

  syms hp
% GenRot=[0;0;1;0;0;0];GenPrism=[0;0;0;0;0;1];
% GenHelical=@(hp)[0;0;1;0;0;hp];
% GetPhiC=@(type)[(type==1)*GenRot+(type==3)*GenCyl];
% GetPhi=@(type,hp)[(type==1)*GenRot+(type==2)*GenPrism+(type==3)*GenHelical(hp)]
% Inertia=@(m,rcm,I)[I+, m*tilde(rcm);tilde(m*rcm),[m 0 0;0 m 0; 0 0 m]];
%% Preparatory steps

for k=1:Nbodies
    rcm(coordRange,k)=str2sym(['[rcm' num2str(k) 'x;rcm' num2str(k) 'y*' num2str(EnableOffsets) ';rcm' num2str(k) 'z*' num2str(EnableOffsets) ';]']);%Center of mass in link frame
    Lp(coordRange,k)=str2sym(['[L' num2str(k) 'x;L' num2str(k) 'y;L' num2str(k) 'z;]']);%Center of mass in link frame
end
for k=1:Nbodies, 
    I0(coordRange,coordRange,k)=GetSimVecSymmetric([3,3],{['I' num2str(k) '_']}); %I(4,4,k)=m(k); I(1:3,4)=tilde(m(k)*rcm(:,k));I(4,1:3)=I(1:3,4) 
    I{k}=rbi(I0(coordRange,coordRange,k),rcm(:,k),m(k));
end

p=[[DH{:,5}]]-1;
jointType=[DH{:,7}];
sigma=jointType;
%%
% % for k=1:Nbodies
% %     % eval(['T(' num2str(DH{k,5}) ',' num2str(DH{k,6}) ').T=DHProximal(DH{k,1},DH{k,2},DH{k,3},DH{k,4})']);
% %     % X(:,:,k,p(k)+1)=BuildX(T(p(k)+1,k).T(1:3,1:3),p);
% %     T{k}=DH2XT(DH{k,1},DH{k,2},DH{k,3},DH{k,4});
% %     [ X{i}, Phi{i} ]=jcalc(jointType(k),q(k));
% %     Phi(:,k)=GetPhi(jointType(k),0);
% % end

%%
% H=zeros(Ndof,Ndof);
pitch(1:Nbodies)=jointType;%Initialize as all rotoidal joint
%% Calculation requires precalculation of transformations between connected bodies
for k=1:Nbodies
%% Initialization of Xtree
 if p(k) == 0
    Xtree{1} = xlt([0 0 0]);
  else
    Xtree{k} = rotx(DH{k,1}) * xlt([-DH{k,2},0,0]);%Initialize link transformation
 end
   % [ XJ, S{k} ] = jcalc(jointType(k), q(k),0);
  XJ = rotz((1-sigma(k))*-q(k)) * xlt([0 0 sigma(k)*q(k)]);
  Phi{k} = str2sym(['[0;0;1;0;0;' num2str(pitch(k)) ']']);  
  % vJ = S{k}*qd(k);%Speed
  Xup{k} = XJ * Xtree{k};%Get joint transformation as composition of a inlink transformation and a cross/link transformation
end
%% JSIM calculation loop (CRBA algorithm)


for i=1:Nbodies
    IC{i}=I{i};%Subtree Inertia 6x6 matrix initialization
end

%% Precalculation 
% for i = Nbodies:-1:1
%   if p(i) ~= 0
%     IC{p(i)} = IC{p(i)} + Xup{i}'*IC{i}*Xup{i};
%   end
% end
%% Loop over tree
% loop over tree bodies
disp(['============>>> SYSTEM DESCRIPTION <<<==================']);
disp(['Number of dof: ' num2str(Ndof)]);
disp(['Number of bodies: ' num2str(Nbodies)]);
disp(['Used DH table: ']);
disp(DH);
disp(['--------------------------------------------------------'])
disp('Loop over tree: calculating JSIM mass matrix using CRBA algorithm....')
tp=tic;
for i=Nbodies:-1:1
    tb=tic;
    disp(['Evaluating body ' num2str(i) ]);
    F=IC{i}*Phi{i};%Initialize utility variable F
    H(i,i)=transpose(Phi{i})*F;%Define diagonal terms
    if p(i)~= 0
        IC{p(i)}= IC{p(i)}+transpose(Xup{i})*IC{i}*Xup{i};
    end
    j=i;
    while p(j)~= 0
        F=simplify(transpose(Xup{j})*F,'Seconds',SimplificationTime);
        j=p(j);
        H(i,j)=simplify(transpose(F)*[Phi{j}],'Seconds',SimplificationTime);
        %H(j,i)=transpose(H(i,j));
        H(j,i)=H(i,j);
    end%while
    tprocb(Nbodies-i+1)=toc(tb);
    plot(Nbodies-[i:Nbodies]+1,log(tprocb),'-bs');grid on;xlabel('Body #');ylabel('log computational time[s]');
    pause(0.001);
end

% % %Use simmetry
% % for i=1:Nbodies
% %     for j=i:Nbodies
% %         H(i,j)=H(j,i);
% %     end
% % end
tproc=toc(tp)
prof=profile('info');
disp(prof)

if SaveResults
disp('Saving results...')
for i=1:Ndof
    for j=1:Ndof
    tmp=H(i,j);
    save(['Session_' RunID '_' num2str(Ndof) '_H_' num2str(i) '_' num2str(j) ],'tmp');
    end
end
% save(['Session_' RunID '_' num2str(Ndof)]);
disp('Done.')
end

%% Post processing and statistics
Vars=symvar(H)
NumOfVars=length(Vars)
% NumOfTerms=length(children(H))
%% Mass matrix simplifications
% profile viewer
% p = profile('info')
% save myprofiledata p
diary off
