%% GDA geometric dynamics algorithm
% Soldà Enrico
%  Last modified 17/08/2024
%  Based on article: "Reducing the computational complexity of inertia-matrix
%  calculation for high DOF robots" from Mohammad Safeea, Richard Bearee
%
% Method is based on explotation of geometrical properties of open chain
% and direct estimation of koint to joint accelerations induces
% Used and useful for benchmarking of different formulations of inertia
% matrix with alterantive methods.
% Geometrc method is considerably faster than analytical methods both for
% jacobians calculation as well for calculation of other quantities and do
% not requires any differentiation as Newton Euler method as well or
% customized Euler Lagrange methods
%
%  This function is used to calculate the joint space inertia matrix (JSIM)
%  of a serial open chain mechanism (from manipulator concept) with only revolute 
%  joints using GDA algorithm.
%  Modified Denavit Hartenberg convention is used

function [At]=GetMassMatrixGDA(T,Pcii,Icii,mcii)
% At: is the joint space inertia matrix
% T: is 4x4xn matrix, representing the homogeneous transformations of the
% link frames
% Thus matrix T(:,:,i) represents 4x4 homogeneous transform of frame i with
% respect to reference frame
% Pcii: is 3xn matrix, each column Pcii(:,i) represent the coordinate
% vector of the center of mass of link i in the local frame of that link.
% Icii: is 3x3xn matrix, thus Icii(:,:,i) matrix represents the
% inertia of link i, around its center of mass represented in frame i
% mcii: is nx1 column vector, while each element mcii(i) represents 
% the mass of link i

% Adaptation to symbolic need to modify variables initialization
n=max(size(mcii)); 
% Ai=zeros(3,n,n);Ci=zeros(3,n,n);
% Ici=zeros(3,3);
% L=zeros(3,3);
% Pcii_A=zeros(3,n);
for i=1:n
    Pcii_A(:,i)=T(1:3,1:3,i)*Pcii(:,i);
    Pci=Pcii_A(:,i)+T(1:3,4,i);
    Ici=T(1:3,1:3,i)*Icii(1:3,1:3,i)*T(1:3,1:3,i)';
    for j=1:i
    Ai(:,j,i)=Ici*T(1:3,3,j);
    Pcij=Pci-T(1:3,4,j);
    Ci(:,j,i)=cross1(T(1:3,3,j),Pcij);
    end
end
% Adaptation to symbolic need to modify variables initialization
% At=zeros(n,n);Bt=zeros(n,1);
% Fac_C=zeros(3,n);
% Mac_A=zeros(3,n);
% Pjp1_j=zeros(3,1);

start=n-1;
j=n;

for k=1:n 
     Mac_A(:,k)=Ai(:,k,j)+cross1(mcii(j)*Pcii_A(:,j),Ci(:,k,j));            
end
Fac_C=mcii(j)*Ci(:,:,j);         
At(j,:)=T(1:3,3,j)'*Mac_A;
    
  
for j=start:-1:1
    Pjp1_j=T(1:3,4,j+1)-T(1:3,4,j);
        for k=1:j 
            Mac_A(:,k)=Mac_A(:,k)+Ai(:,k,j)+cross1(Pjp1_j,Fac_C(:,k))+cross1(mcii(j)*Pcii_A(:,j),Ci(:,k,j));            
        end
    Fac_C(:,1:j)=Fac_C(:,1:j)+mcii(j)*Ci(:,1:j,j);        
    At(j,1:j)=T(1:3,3,j)'*Mac_A(:,1:j);
    
end
for i=1:n
    for j=1:n
        if i>j
            At(j,i)=At(i,j);
        end
    end
end
end
%% Cross product calculation
%%
function c=cross1(a,b)
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:);
     a(3,:).*b(1,:)-a(1,:).*b(3,:);
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
end