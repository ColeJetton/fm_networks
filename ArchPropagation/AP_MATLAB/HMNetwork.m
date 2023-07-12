function [A] = HMNetwork(n,m,r,p,q)

% This function takes produces a matrix A which is a hierarchically modular
% network of n nodes, m modules per level, and r levels where p is
% the edge probability within a module and q is the connectivity off the
% hierarchy
% 
% Keep track of versions here: 
% Date: Version 1: 9 October 2015
% Author: Andy Dong

% Initial lowest level hierarchy module
% Number of nodes in lowest level hierarchy
lnodes=n/(2*power(2,r-1))
% Lowest level hierarchy module
A=modular(lnodes,m,p,lnodes/m,q,1)
for i=1:r,
    p*power(q,i)
    P=RandNetwork(power(2,i-1)*lnodes,p*power(q,i))
    Aprime=[A,P;P,A]
    A=Aprime;
end

A = triu(A);
A = A+A';
A(logical(eye(size(A))))=0;