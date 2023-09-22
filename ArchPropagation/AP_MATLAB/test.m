clc; clear;
format compact
tic %note! Takes roughly 6 minutes to run
% Number of components in the system
nodes=1024;
% Reference system
modules=1;

% Ideal Modular
% b=modular(nodes,modules,1,nodes/modules,0,1);

% Modular
 b=modular(nodes,modules,.5,nodes/modules,.125,1);



% Hierarchically Modular
% Number of levels of hierarchy (should be >= 2)
% levels = 4;
%b=modular(nodes,modules,.9,nodes/modules,.5,.5);

% Antenna case study
% load('antennacpm');
% Convert to a component DSM
% b=directlikelihood;
% b(b>0)=1;
b=triu(b);


b=b+b';

b(logical(eye(size(b))))=0;
figure(2)
imshow(b)