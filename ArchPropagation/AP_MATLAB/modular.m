function [cij]=modular(n,nm,pm,km,pc,pe)
%*****************************************
%
%   This function generates a modular network of size n.
%
%   Inputs: n size of network
%           nm number of modules
%           pm probability of connection in module
%           km size of modules
%           pc probability of intermodule connections
%           pe probability connection exponent
%              pe is used to scale pc
%
%   Outputs: cij modular network
%
%	Author: Richard Gray, School of Physics, The University of Sydney,
%		    N.S.W 2006, Australia
%   Using:  Matlab 7.0

%   Modifications: (Somwrita Sarkar, 2011): Instead of pe being used to
%   scale pc, pc is intermodule connection probablity, and p sets the rate
%   of decay for pc. Further, instead of pc*pm being the intermodule, use
%   pc directly - for cleaner analytics (values of eigenvalues). 
%
%           5/03/2007 
%   
%   Easily generalized to case where inter module connections are created
%   with probability that varies with level.
%
%--------------------------------------------------------------------------

cij=zeros(n,n);

% Do some checking

if (n~=nm*km) disp('Parameters n,m, and km incompatable'); return; end;
%if (kc > km^2) disp('Number interconnections to large'); return; end

levels=log2(nm);
if (levels-round(levels)~=0.0) disp('Number of modules not a power of 2'); return; end; 

% First generate the randomly connected modules. This is done
% by using connectivity.m to make a random network of the correct 
% size and then moving each module to the correct position in cij.

    
for a=1:nm
    module=connectivity(km,'1a',pm);  % create module
    cij(km*(a-1)+1:km*a, km*(a-1)+1:km*a)=module; % put module in correct position

end;

% Now add the connections between modules. Again use connectivity.m
% to generate a random matrix of the correct size and probability of
% connection. Then move to the correct position in cij. As loop over the 
% level number the size of the intermodule matrix increases while
% the number of the connections remains constant. 
% One intermodule connection is added to ensure that the network
% is strongly connected (if pc is too small there is a chance no
% connections are added.)

%impc=pc*pm;
impc = pc; 

for a=1:levels
    
    % Below the diagonal
    
    for b=1:nm/(2^a)
        
        imsize=km*2^(a-1);

        % create intermodule connection matrix
        intermod=connectivity(imsize,'1a',impc);
 
        %if sum(intermod) == 0 
        %    intermod=connectivity(imsize,'5b',1); % make sure there is at least
        %end                                       % one connection      
        
        % move connections into position
        cij(2^(a-1)*km+1+(b-1)*2^a*km:2^a*km+(b-1)*2^a*km,...
        1+(b-1)*2^a*km:2^(a-1)*km+(b-1)*2^a*km)=intermod; 

    end;
    
    % Above the diagonal
    
    for b=1:nm/(2^a)
        
        imsize=km*2^(a-1);
        
        % create intermodule connection matrix
        intermod=connectivity(imsize,'1a',impc);
        
        
        %if sum(intermod) == 0 2^(a-1)*km+1+(b-1)*2^a*km:2^a*km+(b-1)*2^a*km;
        %    intermod=connectivity(imsize,'5b',1);  % make sure there is at least
        %end                                        % one connection        
        
        % move connections into position
        cij( 1+(b-1)*2^a*km:2^(a-1)*km+(b-1)*2^a*km,...
        2^(a-1)*km+1+(b-1)*2^a*km:2^a*km+(b-1)*2^a*km)=intermod;

    end;
    
    impc=impc*pe;
    %impc=impc*(4^(pe-1)); % scale the intermodule connection probability 
                          % for the next level  
end;

% Display cij

%colormap(gray)

%subplot(2,2,1)
%imagesc(cij);
%axis square;
%title('Modular Connectivity Matrix - C(N)');

%subplot(2,2,2)
%[order] = corrmap2(cij);
%imagesc(cij(order,order),[0 max(max(cij-diag(diag(cij))))]);
%axis square;
%title('Reordered Modular Connectivity Matrix - C(N)');
    
