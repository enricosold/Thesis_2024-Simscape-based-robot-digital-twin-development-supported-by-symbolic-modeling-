function [EOM,vars]=Hamilton(H,p,q)
% Soldà Enricò
% Last modified 6 August 2024
% Calculate equations of motion (homogeneous part) for a given Hamiltonian function and given
% generalized momenta and coordinates vectors used in its definition
% Generalized speed and coordinates shall be defined as symbolic function of time
% for proprer functioning
% For complex Hamiltonians, usually when dependance on function of coordinates is included in kinetic terms,
% process may become time consuming
% Generalized variables and speeds shall be in same order

vars=symvar(H);
syms t;
te=tic;
for k=1:length(q)
    EOM(k,1)=-diff(H,q(k))
    EOM(k,2)=diff(H,p(k));
end
tp=toc(te);
fprintf('Done in %f',tp);
end