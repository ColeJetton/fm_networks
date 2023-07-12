function [S] = AntennaPerturbMatrix(A,c,pw,pu)

% This function takes in A and c (component being changed), and produces a matrix S starting from
% a random matrix with probability p and then rewired such that
% S(A==0)=1 with probability pw (probability of wiring)
% S(A==1)=-1 with probability pd (probability of unwiring)
% A is assumed square
% 
% Keep track of versions here: 
% Date: Version 1: 10 November 2015
% Author: Andy Dong
S=zeros(size(A,2));

% Where there is an edge in A with component c perturbation may only be 0 or -1
[~,vu2]=find((A(c,:)+S(c,:))==1);

Dis=(rand(length(vu2),1)>=pu); %pairs to disconnect

for i=1:length(Dis),
    if(Dis(i)==0) % Fails probability test
        S(c,vu2(i))=0; % Make no perturbation
    else
        S(c,vu2(i))=-1; % Perturb by removing edge
    end
end

% Where there is no edge in A with component c perturbation may only be 0 or 1
[~,vw2]=find((A(c,:)-S(c,:))==0);

Con=(rand(length(vw2),1)>=pw);%pairs to connect

for i=1:length(Con),
    if(Con(i)==0) % Fails probability test
        S(c,vw2(i))=0; % Make no perturbation
    else
        S(c,vw2(i))=1; % Perturb by adding edge
    end
end

% Make S symmetric to prevent imaginary eigenvalues
S = triu(S);
S = S+S';
S(logical(eye(size(S))))=0;