function  Tree = autoTree( nb, bf, skew, taper )
% Soldà Enrico
%
% Featherstone R. "Rigid Body Dynamics Algorithms" textbook is referred
%
% For compatibility with other functions we use a structure for tree
% description but alternative representations are used in other places
%
% autoTree  Create System Models of Kinematic Trees
% autoTree(nb,bf,skew,taper)  creates system models of kinematic trees
% having revolute joints, mainly for the purpose of testing dynamics
% functions.  nb and bf specify the number of bodies in the tree, and the
% branching factor, respectively.  The latter is the average number of
% children of a nonterminal node, and must be >=1.  bf=1 produces an
% unbranched tree; bf=2 produces a binary tree; and non-integer values
% produce trees in which the number of children alternates between
% floor(bf) and ceil(bf) in such a way that the average is bf.  Trees are
% constructed (and numbered) breadth-first.  Link i is a thin-walled
% cylindrical tube of   gth l(i), radius l(i)/20, and mass m(i), lying
% between 0 and l(i) on the x axis of its local coordinate system.  The
% values of l(i) and m(i) are determined by the tapering coefficient:
% l(i)=taper^(i-1) and m(i)=taper^(3*(i-1)).  Thus, if taper=1 then
% m(i)=l(i)=1 for all i.  The inboard joint axis of link i lies on the
% local z axis, and its outboard axis passes through the point (l(i),0,0)
% and is rotated about the x axis by an angle of skew radians relative to
% the inboard axis.  If the link has more than one outboard joint then they
% all have the same axis.  If skew=0 then the mechanism is planar.  The
% final one, two or three arguments can be omitted, in which case they
% assume default values of taper=1, skew=0 and bf=1.
% Works symbolically or numerically as well.

if nargin < 4
  taper = 1;
end
if nargin < 3
  skew = 0;
end
if nargin < 2
  bf = 1;
end

if bf<1
    error('Branching factor bf is assumed to be grater or equal 1 (unbranched tree/serial chain)')
end

Tree.NB = nb;
Tree.Ndof=nb;
Tree.Njoints=nb;
len=sym('l%d%d',[nb,3],'real');
rcm=sym('rcm%d%d',[nb,3],'real');
m=sym('m%d',[nb,1],'real');
%% Define joints variables
q=sym('q%d',[nb,1],'real');
dq=sym('dq%d',[nb,1],'real');
ddq=sym('ddq%d',[nb,1],'real');
%% Preliminary setup and initializations
 Tree.geom_prms=[];
 Tree.inertial_prms=[];
 Tree.standard_inertial_prmsvec=[];
%% Tree creation loop

for i = 1:nb
  Tree.jtype{i} = 'R';
  Tree.parent(i) = floor((i-2+ceil(bf))/bf);
  if Tree.parent(i) == 0
    Tree.Xtree{1} = xlt([0 0 0]);
  else
    Tree.Xtree{i} = rotx(skew) * xlt([len(Tree.parent(i),1:3)]);
     Tree.geom_prms=[Tree.geom_prms,len(Tree.parent(i),1:3)];%Patch
  end
  % len(i) = taper^(i-1);
  % mass = taper^(3*(i-1));
  mass=m(i);
  CoM = rcm(i,1:3);
  % Icm = mass * len(i)^2 * diag([0.0025,1.015/12,1.015/12]);
  % Icm = mass * len(i)^2 * diag(sym('I%d',[1,3],'real'));%Custom inertial diagonal model
  Icm = triu(sym(['I' num2str(i) '_%d%d'],[3,3],'real'));%General inertial full tensor model
  Icm=(Icm+Icm')-diag(diag(Icm));%Make it symmetric
  Tree.I{i} = mcI( mass, CoM, Icm );
  Tree.CoM{i}=CoM;
  Tree.m{i}=mass;
  Tree.standard_inertial_prms0{i}=[transpose(diag(Icm)) Icm(1,2) Icm(1,3) Icm(2,3) mass*CoM mass];%Non conventional inertial parameters definition
  Tree.standard_inertial_prmsvec=[Tree.standard_inertial_prmsvec transpose(diag(Icm)) Icm(1,2) Icm(1,3) Icm(2,3) mass*CoM mass];%Non conventional inertial parameters definition
  Tree.standard_inertial_prms{i}=[transpose(diag(Icm)) Icm(1,2) Icm(1,3) Icm(2,3) mass];%Non conventional inertial parameters definition
  Tree.inertial_prms=[ Tree.inertial_prms Tree.standard_inertial_prms{i}];
end
Tree.q=q;
Tree.dq=dq;
Tree.ddq=ddq;
Tree.f_ext=[];
Tree.gravity=str2sym('[g0x;g0y;g0z;dw0x;dw0y;dw0z]');
% drawing instructions

Tree.appearance.base = ...
  { 'box', [-0.2 -0.3 -0.2; 0.2 0.3 -0.06] };

p0 = -1;
for i = 1:nb
  p1 = Tree.parent(i);
  tap = taper^(i-1);
  if p1 == 0
    ptap = 1;
  else
    ptap = taper^(p1-1);
  end
  if ( p1 > p0 )
    Tree.appearance.body{i} = ...
      { 'cyl', [0 0 0; 1 0 0]*tap, 0.05*tap, ...
        'cyl', [0 0 -0.07; 0 0 0.07]*ptap, 0.08*ptap };
    p0 = p1;
  else
    Tree.appearance.body{i} = ...
      { 'cyl', [0 0 0; 1 0 0]*tap, 0.05*tap };
  end
end
Tree.geom_prms=transpose(Tree.geom_prms);
Tree.AdjMatrix=full(sparse([Tree.parent(2:nb),2:nb],[2:nb,Tree.parent(2:nb)],ones(2*length(Tree.parent(2:nb)),1)));
Tree.Graph=graph(Tree.AdjMatrix);