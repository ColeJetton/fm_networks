function [S] = HMPerturbMatrix(A,n,m,p,q,pw,pu)

% This function takes in HM matrix A and produces a matrix S starting from
% a modular matrix with intra-module connection probability p and inter-module
% connection probability q then rewired such that
% S(A==0)=1 with probability pw (probability of wiring)
% S(A==1)=-1 with probability pd (probability of unwiring)
% A is assumed square and the parameters n,m,r should be correct for A

% This is used to perturb a HM network within or outside a module.
% To perturb intra-module set q = 0
% To perturb inter-module set p = 0
% 
% Keep track of versions here: 
% Date: Version 1: 24 November 2015

% Initial lowest level hierarchy module
% Number of nodes in lowest level hierarchy
% lnodes=n/(2*power(2,r-1));
% Lowest level hierarchy module
% S=modular(lnodes,m,p,lnodes/m,q,1);
% for i=1:r,
%     P=RandNetwork(power(2,i-1)*lnodes,p*power(q,i));
%     Sprime=[S,P;P,S];
%     S=Sprime;
% end

S=modular(n,m,p,n/m,q,.5);

% Where there is an edge in A perturbation may only be 0 or -1
[vu1,vu2]=find(A+S==2);
Dis=(rand(length(vu1),1)>=pu); %pairs to disconnect
for i=1:length(Dis),
    if(Dis(i)==0) % Fails probability test
        S(vu1(i),vu2(i))=0; % Make no perturbation
    else
        S(vu1(i),vu2(i))=-1; % Perturb by removing edge
    end
end

% Where there is no edge in A perturbation may only be 0 or 1
[vw1,vw2]=find(A-S==-1);
Con=(rand(length(vw1),1)>=pw);%pairs to connect
for i=1:length(Con),
    if(Con(i)==0) % Fails probability test
        S(vw1(i),vw2(i))=0; % Make no perturbation
    else
        S(vw1(i),vw2(i))=1; % Perturb by adding edge
    end
end

S = triu(S);
S = S+S';
S(logical(eye(size(S))))=0;