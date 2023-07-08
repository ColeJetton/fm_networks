function [A] = RandNetwork(n,p)

A = rand(n,n) < p;
A = triu(A,1);
A = A+A';
% for i = 1:n
%     for j = 1:n
%         if A1(i,j) == 1
%             A(i,j) = 1;
%         end
%     end
% end