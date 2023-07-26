% quick computational experiment to see if we can reduce the time 

% Number of eigenvectors to retain
k=16;
% Number of computation loops
g=1000;
% Memory initialization
a1=cell(1,g);
x1=cell(1,g);
y1=cell(1,g);
z1=cell(1,g);
% Antenna case study
load('antennacpm'); 
% Convert to a component DSM
b=directlikelihood;
b(b>0)=1;
b=triu(b);
b=b+b';
b(logical(eye(size(b))))=0;

for i=1:g
        Er1=AntennaPerturbMatrix(b,3,.1,.1);
        Er2=AntennaPerturbMatrix(b,3,.5,.5);
        Er3=AntennaPerturbMatrix(b,3,.9,.1);
        Er4=AntennaPerturbMatrix(b,3,.1,.9);
        a1{i}=ComputeSine(b,Er1,k)';
        x1{i}=ComputeSine(b,Er2,k)';
        y1{i}=ComputeSine(b,Er3,k)';
        z1{i}=ComputeSine(b,Er4,k)';
end
a1=cell2mat(a1);
x1=cell2mat(x1);
y1=cell2mat(y1);
z1=cell2mat(z1);
%% other sample size 100
% Number of computation loops
g=100;
% Memory initialization
a2=cell(1,g);
x2=cell(1,g);
y2=cell(1,g);
z2=cell(1,g);

for i=1:g
        Er1=AntennaPerturbMatrix(b,3,.1,.1);
        Er2=AntennaPerturbMatrix(b,3,.5,.5);
        Er3=AntennaPerturbMatrix(b,3,.9,.1);
        Er4=AntennaPerturbMatrix(b,3,.1,.9);
        a2{i}=ComputeSine(b,Er1,k)';
        x2{i}=ComputeSine(b,Er2,k)';
        y2{i}=ComputeSine(b,Er3,k)';
        z2{i}=ComputeSine(b,Er4,k)';
end

a2=cell2mat(a2);
x2=cell2mat(x2);
y2=cell2mat(y2);
z2=cell2mat(z2);
%% other sample size 50
% Number of computation loops
g=50;
% Memory initialization
a3=cell(1,g);
x3=cell(1,g);
y3=cell(1,g);
z3=cell(1,g);

for i=1:g
        Er1=AntennaPerturbMatrix(b,3,.1,.1);
        Er2=AntennaPerturbMatrix(b,3,.5,.5);
        Er3=AntennaPerturbMatrix(b,3,.9,.1);
        Er4=AntennaPerturbMatrix(b,3,.1,.9);
        a3{i}=ComputeSine(b,Er1,k)';
        x3{i}=ComputeSine(b,Er2,k)';
        y3{i}=ComputeSine(b,Er3,k)';
        z3{i}=ComputeSine(b,Er4,k)';
end

a3=cell2mat(a3);
x3=cell2mat(x3);
y3=cell2mat(y3);
z3=cell2mat(z3);

%% Statistical comparison

H = zeros(16,12);
P = H;
for i = 1:16
    [H(i,1),P(i,1)]=ttest2(a1(i,:),a2(i,:));
    [H(i,2),P(i,2)]=ttest2(a1(i,:),a3(i,:));
    [H(i,3),P(i,3)]=ttest2(a2(i,:),a3(i,:));
    [H(i,4),P(i,4)]=ttest2(x1(i,:),x2(i,:));
    [H(i,5),P(i,5)]=ttest2(x1(i,:),x3(i,:));
    [H(i,6),P(i,6)]=ttest2(x2(i,:),x3(i,:));
    [H(i,7),P(i,7)]=ttest2(y1(i,:),y2(i,:));
    [H(i,8),P(i,8)]=ttest2(y1(i,:),y3(i,:));
    [H(i,9),P(i,9)]=ttest2(y2(i,:),y3(i,:));
    [H(i,10),P(i,10)]=ttest2(z1(i,:),z2(i,:));
    [H(i,11),P(i,11)]=ttest2(z1(i,:),z3(i,:));
    [H(i,12),P(i,12)]=ttest2(z2(i,:),z3(i,:));
end
mean_P = mean(mean(P));
mean_H = mean(mean(H));

%% stats two

H = zeros(16,4);
P = H;
for i = 1:16
    [H(i,1),P(i,1)]=ttest2(a2(i,:),a3(i,:));
    [H(i,2),P(i,2)]=ttest2(x2(i,:),x3(i,:));
    [H(i,3),P(i,3)]=ttest2(y2(i,:),y3(i,:));
    [H(i,4),P(i,4)]=ttest2(z2(i,:),z3(i,:));
end
mean_P = mean(mean(P))
mean_H = mean(mean(H))