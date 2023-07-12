function [G] = idealModular(n, m)

% Function idealModular generates perfect modular networks with fully
% connected equally sized modules. 
% n: number of nodes in the network
% m: number of modules #g
% size of module: n/m needs to be an integer

g = cell(1,m);
[g{:}] = deal(ones(n/m,n/m));
G = blkdiag(g{:});
G(logical(eye(size(G))))=0;