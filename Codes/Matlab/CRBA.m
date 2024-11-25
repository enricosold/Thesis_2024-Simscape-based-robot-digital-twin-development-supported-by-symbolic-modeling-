function  [H] = CRBA( model)
% Soldà Enrico
% Last modified 13/10/2024
% Featherstone R. "Rigid Body Dynamics Algorithms" textbook is referred
%
% JS_IM Calculate joint space inertia matrix JSIM of equation of motion
% [H]=JS_IM(model) using CRBA(Composite-Rigid-Body) for H.


%% Preliminary system kinematics definition
q=model.q;
for i = 1:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
end
%% CRBA algorithm for trees, calculation of inertia matrix
IC = model.I;				% composite inertia calculation initialization
% Backward passage from tips to root
for i = model.NB:-1:1
    if model.parent(i) ~= 0
        IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i};
    end
end

%H = zeros(model.NB);

for i = 1:model.NB
    fh = IC{i} * S{i};
    H(i,i) = S{i}' * fh;
    j = i;
    while model.parent(j) > 0
        fh = Xup{j}' * fh;
        j = model.parent(j);
        H(i,j) = S{j}' * fh;
        H(j,i) = H(i,j);
    end
end
end