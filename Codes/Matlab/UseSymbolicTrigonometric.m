function [Tmp,vars]= UseSymbolicTrigonometric(M,Ndof)
% Soldà Enrico
% Last modfied  13 Novemebr 2024
% Last modified 18 August 2024
% Vecorially redefine trigonometric functions inside expression (assumed to
% depend only on generalized coordiantes q)
% Tought to be used for inertia matrix but could be use (wary of) for
% general expressions
% Assumed to depend on Ndof coordinates q in number of number of rows of M
if nargin<2
Ndof=size(M,1);
end

q=sym('q%d',[Ndof,1],'real');%Redefine variables vector
c=sym('c%d',[Ndof,1],'real');%Define symbolic cosine
s=sym('s%d',[Ndof,1],'real');%Define symbolic sine
cc=sym('c%d%d',[Ndof,Ndof],'real');%Define symbolic cos of sum of pairs of variables
ss=sym('s%d%d',[Ndof,Ndof],'real');%Define symbolic cos of sum of pairs of variables
vars=[c;s;[cc(:)];[ss(:)]];
Tmp=M;%Initialize recursive substitution with original matrix

for k=1:Ndof
    Tmp=subs(Tmp,cos(q(k)),c(k));
    Tmp=subs(Tmp,sin(q(k)),s(k));
end
% Pairs
for k=1:Ndof
    for h=1:Ndof
        Tmp=subs(Tmp,cos(q(k)+q(h)),cc(k,h));
        Tmp=subs(Tmp,sin(q(k)+q(h)),ss(k,h));
    end
end
% Multiple
        for h=1:Ndof
            for k=h:Ndof
            c_n=sym(['c' strrep(num2str(h:k),' ','')]);
            s_n=sym(['s' strrep(num2str(h:k),' ','')]);
            Tmp=subs(Tmp,cos(sum(q(h:k))),c_n(1));
            Tmp=subs(Tmp,sin(sum(q(h:k))),s_n(1));
            end
        end
    

end