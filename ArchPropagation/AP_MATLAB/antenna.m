clear;
% Number of eigenvectors to retain
k=16;
% Number of computation loops
g=1000;
% Memory initialization
a=cell(1,g);
x=cell(1,g);
y=cell(1,g);
z=cell(1,g);
% Antenna case study
load('antennacpm'); 
% Convert to a component DSM
b=directlikelihood;
b(b>0)=1;
b=triu(b);
b=b+b';
b(logical(eye(size(b))))=0;
%%
% Peturbation Matrix
% Engines is component 3
% Engine Auxiliaries is component 10
% Auxiliary Electrics is component 13
% Cabling & Piping is component 14
Er1=AntennaPerturbMatrix(b,3,.1,.1);
Er2=AntennaPerturbMatrix(b,3,.5,.5);
Er3=AntennaPerturbMatrix(b,3,.9,.1);
Er4=AntennaPerturbMatrix(b,3,.1,.9);
for i=1:g,
        Er1=AntennaPerturbMatrix(b,3,.1,.1);
        Er2=AntennaPerturbMatrix(b,3,.5,.5);
        Er3=AntennaPerturbMatrix(b,3,.9,.1);
        Er4=AntennaPerturbMatrix(b,3,.1,.9);
        a{i}=ComputeSine(b,Er1,k);
        x{i}=ComputeSine(b,Er2,k);
        y{i}=ComputeSine(b,Er3,k);
        z{i}=ComputeSine(b,Er4,k);
end


%Find the mean perturbations
% This calculates the mean perturbation per value of p
for i=1:1,
am(i,:)=mean(cell2mat(a(i,:)'),1);
xm(i,:)=mean(cell2mat(x(i,:)'),1);
ym(i,:)=mean(cell2mat(y(i,:)'),1);
zm(i,:)=mean(cell2mat(z(i,:)'),1);
end

agrandmean=mean(am,1);
xgrandmean=mean(xm,1);
ygrandmean=mean(ym,1);
zgrandmean=mean(zm,1);

% Plot the results
h1=figure('Name','Eigenvector Rotation of Antenna System');
subplot(2,2,1);
scatter(1:size(agrandmean,2),agrandmean,'Marker','o','MarkerEdgeColor','red');
xlabel('Eigenvector index');
ylabel('Angle of rotation');
title('Perturbed pw=.1 pu=.1');
subplot(2,2,2);
scatter(1:size(xgrandmean,2),xgrandmean,'Marker','o','MarkerEdgeColor','red');
xlabel('Eigenvector index');
ylabel('Angle of rotation');
title('Perturbed pw=.5 pu=.5');
subplot(2,2,3);
scatter(1:size(ygrandmean,2),ygrandmean,'Marker','o','MarkerEdgeColor','red');
xlabel('Eigenvector index');
ylabel('Angle of rotation');
title('Perturbed pw=.9 pu=.1');% Plot the results
subplot(2,2,4);
scatter(1:size(zgrandmean,2),zgrandmean,'Marker','o','MarkerEdgeColor','red');
xlabel('Eigenvector index');
ylabel('Angle of rotation');
title('Perturbed pw=.1 pu=.9');
% 
% % Cluster results to find community eigenspace and bulk eigenspace
% [idx,C]=kmeans(transpose(agrandmean),2,'Replicates',10,'Distance','cityblock');	
% % Mean angle community eigenspace
% mace=mean(agrandmean(idx==1));  
% % Mean angle bulk eigenspace
% mabe=mean(agrandmean(idx==2));
% % Number of eigenvalues in community eigenspace
% nece=length(find(idx==1));
% % Number of eigenvalues in bulk eigenspace
% nebe=length(find(idx==2));