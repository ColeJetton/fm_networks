% Trial Run 1

n = 1024;
m = 8;
s = n/m;
p = 0.9;
q = 0.1;
A = modular(n,m,p,s,q,1);
A = triu(A);
A = A+A';
% for i = 1:n
%     A(i,i) = 0;
% end
A(logical(eye(size(A))))=0;
A = sparse(A);
A = A*2;

% Generate Error matrix, a random bernoulli matrix with 0.5 probability of +1 and -1.
% Random network with p=0.5
H = RandNetwork(n,0.5);

% Initialize Error matrix Er
Er = ones(n,n);

% Find all zeros and fill with -1s
% [I,J] = find(H==0);
% for i = 1:length(I)
% Er(I(i),J(i))=-1;
% end
Er(H==0)=-1;

% Set diagonal to zero
%Er = Er - diag(diag(Er));
Er(logical(eye(size(Er))))=0;

Er = sparse(Er);

for i = 1:50
    S(i,:) = ComputeSine(A,Er);
end

m = mean(S);
figure;
scatter(1:100, m);