function [S] = ComputeSineReverse(A,E,k)

% This function takes in A and E, and produces a vector S which is the
% angles by which eigenvectors of A rotate under perturbation by E
% with k large eigenvectors to retain
% Keep track of versions here: 
% Date: Version 1: 1 Sept 2015
% Author: Somwrita Sarkar 
tic
% Compute eigenvectors for A
[v,~] = eigs(A,k,'sm');

% Compute eigenvectors for A+E. Here we assume that E is already
% appropriately scaled.
%a = A + E;
[v1,~] = eigs(A+E,k,'sm');

% Now compute sine between eigenvectors of A and eigenvectors of A+E
for i = 1:k
    c(i) = dot(v(:,i),v1(:,i))/(norm(v(:,i))*norm(v1(:,i))); % Measure cosine
    s(i) = sqrt(1-(c(i)^2)); % Compute sine
    S(i) = asind(s(i)); % Find angle
end
toc