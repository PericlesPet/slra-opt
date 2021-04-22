clear
close all 
%% generate random test problem
close all 
A = randn(1000,400);
[U,S,V] = svd(A,'econ');
% modify to have (slightly) decaying spectrum 
s0 = diag(S);
% s1: simply add C to first k entries
k = 50; 
s1 = s0; 
s1(1:k) = s1(1:k)+100;
A1 = U*diag(s1)*V';
%{ 
s2: transform constant distribution to exponential 
PDF(index) = initArea * lambda * exp(-lambda * index) 
CDF(index) = initArea * [1 - exp(-lambda * index) ] 
%}
s2 = s0; 
lambda = 0.01;                      % decay rate
alpha = sum(s0) / 50;               % constant factor  
s2(:) = lambda * alpha *  s2(:) .* exp(-lambda * (1:length(s2)))'; 
A2 = U*diag(s2)*V';

subplot(2,3,1)
plot(s0)
title("distribution of s0")
subplot(2,3,2)
plot(s1)
title("distribution of s1")
subplot(2,3,3)
plot(s2)
title("distribution of s2")

subplot(2,3,4)
plot(cumsum(s0))
title("cumulative of s0")
subplot(2,3,5)
plot(cumsum(s1))
title("cumulative of s1")
subplot(2,3,6)
plot(cumsum(s2))
title("cumulative of s2")

%% reconstruct A, A1, A2
[U0,D0,V0] = svd(A);
[U1,D1,V1] = svd(A1);
[U2,D2,V2] = svd(A2);

%% Fixed-rank SVD approx. (known k)
E = U(:,1:k);
F = D(1:k,1:k)* V(:,1:k)';

%% Fixed precision SVD approx. (known epsilon, not known k)
kArray = [];
maxEpsilon = norm(A);                           % Maximum epsilon value practically the norm of A
minEpsilon = maxEpsilon / 100;                  % Minimum Epsilon and epsilon increments are 1% of max

for epsilon = minEpsilon:minEpsilon:maxEpsilon
    k = sum(diag(D) > epsilon);
    kArray = [kArray k];
end

E = U(:,1:k);
F = D(1:k,1:k)* V(:,1:k)';

epsilon = minEpsilon:minEpsilon:maxEpsilon;
plot(kArray, epsilon / norm(A))
ylabel("epsilon")
xlabel("rank k")
%% GRAM SCHMIDT

[m, n] = size(A); 
Q = zeros(m, n); 
R = zeros(n, n);

epsilon = 0.01;
epsilons = []
for j = 1:n                          
    v = A(:,j);                     % v is column j of A 
    for i=1:j-1                      
        R(i,j) = Q(:,i)' * A(:,j);  % modify A(:,j) to v 
        v = v - R(i,j)*Q(:,i);      % subtract the projection (q_i' a_j) * q_i = (q_i' v) * q_i
    end                             % v is perpendicular to q_1...q_(j-1)
    R(j,j) = norm(v);
    Q(:,j) = v / R(j,j);
    epsilons = [epsilons norm(A-Q*Q'*A)];
    if norm(A-Q*Q'*A) < epsilon
        k = j;
        break
    end
end

norm(A-Q*Q'*A)


aaa = sort(diff(epsilons))
