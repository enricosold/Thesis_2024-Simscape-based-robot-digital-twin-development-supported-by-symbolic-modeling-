function  tau = ID0( model)
% Soldà Enrico
%
%
% ID0  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID0(model) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

%% Auxiliary functions definition
% crf  spatial/planar cross-product operator (force).
% crf(v)  calculates the 6x6 (or 3x3) matrix such that the expression
% crf(v)*f is the cross product of the motion vector v with the force
% vector f.  If length(v)==6 then it is taken to be a spatial vector, and
% the return value is a 6x6 matrix.  Otherwise, v is taken to be a planar
% vector, and the return value is 3x3.
crf=@( v )[-crm(v)'];

a_grav = model.gravity;
q=model.q;
qd=model.dq;
qdd=model.ddq;
f_ext=model.f_ext;

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
  f{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
end

if nargin == 5
  f = apply_external_forces( model.parent, Xup, f, f_ext );
end

for i = model.NB:-1:1
  tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i};
  end
end
end
%% Auxiliary functions definition
function  vcross = crm( v )

% crm  spatial/planar cross-product operator (motion).
% crm(v)  calculates the 6x6 (or 3x3) matrix such that the expression
% crm(v)*m is the cross product of the motion vectors v and m.  If
% length(v)==6 then it is taken to be a spatial vector, and the return
% value is a 6x6 matrix.  Otherwise, v is taken to be a planar vector, and
% the return value is 3x3.

if length(v) == 6

  vcross = [  0    -v(3)  v(2)   0     0     0    ;
	      v(3)  0    -v(1)   0     0     0    ;
	     -v(2)  v(1)  0      0     0     0    ;
	      0    -v(6)  v(5)   0    -v(3)  v(2) ;
	      v(6)  0    -v(4)   v(3)  0    -v(1) ;
	     -v(5)  v(4)  0     -v(2)  v(1)  0 ];

else

  vcross = [  0     0     0    ;
	      v(3)  0    -v(1) ;
	     -v(2)  v(1)  0 ];
end
end
