function [cij] = connectivity(n,type,param1,param2,param3,param4,param5,param6)
%***********************************************************************

% This function is used to generate a connectivity matrix
% for a network of brain components.

%	Author: Richard Gray, School of Physics, The University of Sydney,
%		    N.S.W 2006, Australia
%   Using:  Matlab 7.0

% Returns a connectivity matrix of size n.
% The input type is the type of connectivity used.
% param1 and param2 are parameters dependent on type
%
%	type
%	----
%
% 	1a	Random
%		param1 = Probability of connection
%
%	1b  Random symmetric (Network is symmetric)
%		param1 = Probability of connection
%
%	1c	Random partially symmetric
%		param1 = Probability of connection
%		param2 = Probability of connection being symmetric
%
%	2a	Small World (Newman/Watts/Monasson  no rewiring)
%		param1 = degree of underlying regular network
%		param2 = probability of shortcut
%
%	2b	Small World (Watts/Strogatz  rewiring)
%		param1 = degree of underlying regular network
%		param2 = probability of rewiring
%
%      2c  Small World (ours rewiring)
%		param1 = degree of underlying regular network
%		param2 = probability of rewiring
%
%	3	Scale Free, future work
%
%	4a	Sigmoid distribution
%       param1 = minimum probability
%       param2 = steepness of transition
%       param3 = position where transition occurs
%
%	5a  Randomly distributed with fixed number of connections
%       No self conections. Used by Sporns et al
%       param1= Number of connections
%
%   5b  Same as 5a but self connections allowed.
%       Used by modular.m
%		param1= Number of connections
%
%   5c  Randomly distributed with fixed number of connections but
%       each node has the same indegree. No self conections. Used by Sporns et al
%       param1= Number of connections
%
%   6a  Modular network. Using modular2.
%       param1 = number of modules
%       param2 = probability of connection in module
%       param3 = size of modules
%       param4 = number of intermodule connections
%       param5 = probability connection exponent
%                pe is used to scale pc
%
%   6b  Partially Symmetric Modular network. Using modular2symm.
%       param1 = number of modules
%       param2 = probability of connection in module
%       param3 = size of modules
%       param4 = number of intermodule connections
%       param5 = probability connection exponent
%                pe is used to scale pc
%       param6 = probability of symmetric connection
%
%-----------------------------------------------------------------------

%intialize
cij=zeros(n,n);

% Generate connectivity matrix

switch type
    case{'1a'}
        temp = rand(n,n);
        [i,j,s] = find(temp <= param1);
        weights = ones(size(s));
        cij = full(sparse(i,j,weights,n,n));
        %cij = double(temp<=param1); %better way to code this

    case{'1b'}
        temp = triu(rand(n,n));
        temp = tril(temp',-1)+temp;
        [i,j,s] = find(temp <= param1);
        weights = ones(size(s));
        cij = full(sparse(i,j,weights,n,n));
    case{'1c'}
        cij = random_near_symmetric(n,param1,param2);
    case{'2a'}
        cij = small_world(n,param1,param2,'1');
    case{'2b'}
        cij = small_world(n,param1,param2,'2');
    case{'2c'}
        cij = small_world(n,param1,param2,'3');
    case {'3'}
        disp('Not ready yet.')
    case{'4a'}
        % Sigmoid distribution
        for i=1:n
            for j=1:n
                d=distance(i,j,n);
                p = (1-param1)/(1+exp(param2*(d-param3))) + param1;
                cij(i,j) = ranbin(p);
            end;
        end;
    case{'5a'}
        % All param1 connections are randomly distributed, diagonal is
        % blocked.
        temp = ones(n)-eye(n);
        [i,j,s] = find(temp);
        r = randperm(n^2-n);
        weights = ones(1,param1);
        cij = sparse(i(r(1,1:param1)),j(r(1,1:param1)),weights,n,n);
    case{'5b'}
        % All param1 connections are randomly distributed, diagonal is NOT
        % blocked.
        temp = ones(n);
        [i,j,s] = find(temp);
        r = randperm(n^2);
        weights = ones(1,param1);
        cij = sparse(i(r(1,1:param1)),j(r(1,1:param1)),weights,n,n);
    case{'5c'}
        % All param1 connections randomly distributed, but we
        % fix in-degree to param1/n
        cols = ones(1,n).*param1/n;
        temp = ones(n)-eye(n);
        cij = zeros(n);
        for unit=1:n
            [ic jc sc] = find(temp(unit,:) == 1);
            rpmake = randperm(length(ic));
            cij(jc(rpmake(1:cols(unit))),unit) = ones(cols(unit),1);
        end;
        cij = sparse(cij);
    case{'6a'}
        cij = modular2(n,param1,param2,param3,param4,param5);
    case{'6b'}
        cij = modular2symm2(n,param1,param2,param3,param4,param5,param6);
    otherwise
        disp('Type of connectivity not specified.')
end

