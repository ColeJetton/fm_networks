function S = RotationPairs(A)
% This function takes in an eigenvalue rotation vector A
% and returns the index value corresponding to the last eigenvalue
% in the community eigenspace
% See 	arXiv:1510.07064 [physics.soc-ph]
% 
% Keep track of versions here: 
% Date: Version 1: 24 November 2015
% Author: Andy Dong
S=[];

% Has the peak been reached?
peak=0;
% Has a trough been reached?
trough=0;
% Allow for some noise in eigenvalues
noise=rad2deg(pi/256);
% Unusual situation for ideal modular network
if A(1) > rad2deg(pi/16),
    peak=1;
end
for i=1:length(A)-1,
    if not(peak),
        if A(i+1)<A(i),
            if A(i+1)+noise<A(i),
                peak=1;
            end
        end
    else
        if A(i+1)>A(i),
            if A(i+1)-noise>A(i) && A(i+1)-noise>A(1),
                trough=1;
            end
        end
    end
    if trough,
        S=i;
        break;
    end
end

if isempty(S), % Oops the eigenvalues increase monotonically into the bulk
    diff=zeros(length(A)-1);
    for i=1:length(A)-1,S
        diff(i)=A(i+1)-A(i);
    end
    % Find the largest jump to the bulk eigenspace
    [~,S]=max(diff(2:end));
    S=S+1;
end