clc; clear;
format compact
tic %note! Takes roughly 6 minutes to run
% Number of components in the system
nodes=1024;
% Erdos-Renyi connection probability
% p=[.9 .8 .7 .6 .5 .4 .3 .2 .1];
% Modular Perturbation Probability
p = [0.05 0.025 0.01 0.005 0.005 0.001];%p = [0.05 0.025 0.01 0.005 0.025 0.001];
% Number of eigenvectors to retain
k=100;
% Number of computation loops
g=100;
% Memory initialization
Pt=cell(1,length(p));
a=cell(length(p),g);
x=cell(length(p),g);
y=cell(length(p),g);
am=zeros(length(p),k);
xm=zeros(length(p),k);
ym=zeros(length(p),k);
% Reference system
modules=8;

% Ideal Modular
% b=modular(nodes,modules,1,nodes/modules,0,1);

% Modular
% b=modular(nodes,modules,.9,nodes/modules,.6,1);


% Hierarchically Modular
% Number of levels of hierarchy (should be >= 2)
% levels = 4;
b=modular(nodes,modules,.9,nodes/modules,.5,.5);

% Antenna case study
% load('antennacpm');
% Convert to a component DSM
% b=directlikelihood;
% b(b>0)=1;
b=triu(b);


b=b+b';

b(logical(eye(size(b))))=0;

%%
for i=1:length(p),
    % Compute perturbation
    for j=1:g,
        % Random Peturbation Matrix
%         Er1=PerturbMatrix(b,p(i),.5,.5);
%         Er2=PerturbMatrix(b,p(i),.9,.1);
%         Er3=PerturbMatrix(b,p(i),.1,.9);
%         Pt{i}=Er1;
        %            Modular Perturbation Matrix
        %         Er1=ModularPerturbMatrix(b,modules,0,p(i),.5,.5);
        %         Er2=ModularPerturbMatrix(b,modules,0,p(i),.9,.1);
        %         Er3=ModularPerturbMatrix(b,modules,0,p(i),.1,.9);
        %         Pt{i}=Er1;
        % Hierarchically Modular Perturbation Matrix
        Er1=HMPerturbMatrix(b,nodes,modules,0,p(i),.2,.8);
        Er2=HMPerturbMatrix(b,nodes,modules,0,p(i),.3,.7);
        Er3=HMPerturbMatrix(b,nodes,modules,0,p(i),.4,.6);
        a{i,j}=ComputeSine(b,Er1,k);
        x{i,j}=ComputeSine(b,Er2,k);
        y{i,j}=ComputeSine(b,Er3,k);
    end
    % Erdos-Renyi random perturbation
    %r=erdosRenyi(nodes,p(i),4);
    %Er=full(r.Adj);
    %Er=triu(Er);
    %Er=Er+Er';
    %Er(logical(eye(size(Er))))=0;
    %Pt{i}=Er;
    % Compute perturbation
    %a(i,:)=ComputeSine(b,Er,k);
    % Somwrita's random perturbation
    %      H=RandNetwork(nodes,p(i));
    %      Er=ones(nodes,nodes);
    %      Er(H==0)=-1;
    %      Er(logical(eye(size(Er))))=0;
    %      Er = sparse(Er);
    %      Pt{i}=Er;
end
% Plot the results
h1=figure('Name','Eigenvector Rotation');
subplot(2,2,1);
%scatter(1:size(brot,2),sort(real(brot),'descend'),'*');
imshow(b);
%xlabel('Eigenvector index');
%ylabel('Angle of rotation');
title('Unperturbed System');
%subplot(2,2,2);
%imshow(b+Pt{1});
%title('Perturbed System p=.9');
%markers = ['o','+','.','x','s','d','^','p','h'];
%h2=figure('Name','Perturbed System');
%hold on;
%for i=1:length(p),
%    scatter(1:size(a(i,:),2),sort(real(a(i,:)),'descend'),markers(i));
%end

%Find the mean perturbations
%am=mean(a,1);
% This calculates the mean perturbation per value of p
for i=1:length(p),
    am(i,:)=mean(cell2mat(a(i,:)'),1);
    xm(i,:)=mean(cell2mat(x(i,:)'),1);
    ym(i,:)=mean(cell2mat(y(i,:)'),1);
end

agrandmean=mean(am,1);
xgrandmean=mean(xm,1);
ygrandmean=mean(ym,1);

subplot(2,2,2);
%scatter(1:size(am,2),sort(real(am),'descend'),'Marker','o','MarkerEdgeColor','red');
scatter(1:size(agrandmean,2),agrandmean,'Marker','o','MarkerEdgeColor','red');
%scatter(1:size(a(i,:),2),sort(real(a(i,:)),'descend'),'o');
xlabel('Eigenvector index');
ylabel('Angle of rotation');
%legend('p=.9','p=.8','p=.7','p=.6','p=.5','p=.4','p=.3','p=.2','p=.1');
title('Perturbed pw=.5 pu=.5');
subplot(2,2,3);
scatter(1:size(xgrandmean,2),xgrandmean,'Marker','o','MarkerEdgeColor','red');
xlabel('Eigenvector index');
ylabel('Angle of rotation');
title('Perturbed pw=.9 pu=.1');
subplot(2,2,4);
scatter(1:size(ygrandmean,2),ygrandmean,'Marker','o','MarkerEdgeColor','red');
xlabel('Eigenvector index');
ylabel('Angle of rotation');
title('Perturbed pw=.1 pu=.9');



toc
% Print Results
str=sprintf('%u & %.1f & %.1f & %.1f & %u & %.1f & %.1f & %.1f & %u & %.1f & %.1f & %.1f & %.1f',nace, ...
    mace, mabe, mabe-mace, nxce, mxce, mxbe, mxbe-mxce, nyce, myce, mybe, mybe-myce, max([mabe-mace,mxbe-mxce,mybe-myce])-min([mabe-mace,mxbe-mxce,mybe-myce]));
fprintf(str);