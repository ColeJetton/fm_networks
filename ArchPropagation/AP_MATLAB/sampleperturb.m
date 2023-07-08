clear;
tic
% Number of components in the system
nodes=16384;
% Random connection probability
p=.5;
% Number of eigenvectors to retain
k=100;
% Number of calculations
loops=100;
% Reference system
modules=16;
% Memory allocation
a1=zeros(loops,k);
a2=zeros(loops,k);
a3=zeros(loops,k);
% Modular
b1=modular(nodes,modules,.1,nodes/modules,.1,1);
b1=triu(b1);
b1=b1+b1';
b1(logical(eye(size(b1))))=0;
b2=modular(nodes,modules,.9,nodes/modules,.6,1);
b2=triu(b2);
b2=b2+b2';
b2(logical(eye(size(b2))))=0;
% Random network perturbation
for i=1:loops,
    Er1=PerturbMatrix(b1,.5,.5,.5);
    Er2=PerturbMatrix(b2,.5,.5,.5);
    a1(i,:)=ComputeSine(b1,Er1,k);
    a2(i,:)=ComputeSine(b2,Er2,k);
    a3(i,:)=ComputeSineReverse(b1,Er1,k);
end

% Plot the results
h1=figure('Name','Eigenvector Rotation of Sample Networks');
subplot(1,2,1);
scatter(1:size(a1,2),mean(a1),'Marker','o','MarkerEdgeColor','red');
xlabel('Eigenvector index');
ylabel('Angle of rotation');
title('Sparse Modular Network');
subplot(1,2,2);
scatter(1:size(a2,2),mean(a2),'Marker','o','MarkerEdgeColor','red');
xlabel('Eigenvector index');
ylabel('Angle of rotation');title('Perturbation of Sparse Modular Network');
title('Dense Modular Network');

% Cluster results to find community eigenspace and bulk eigenspace
[idx,C]=kmeans(transpose(mean(a)),2,'Replicates',10,'Distance','cityblock');	
% Mean angle community eigenspace
mace=mean(a(idx==1));  
% Mean angle bulk eigenspace
mabe=mean(a(idx==2));
% Number of eigenvalues in community eigenspace
nece=length(find(idx==1));
% Number of eigenvalues in bulk eigenspace
nebe=length(find(idx==2));
toc