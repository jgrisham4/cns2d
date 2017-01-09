%% Clearing workspace

clc,clear, close all

%% Inputs

d  = [ 5 4 2 1 8 ];
ld = [1 5 -10 2];

ud = [-2 5 -10 1];
sol_act = ones(numel(d),1);

% Forming A-matrix
A = diag(d);
A = A + diag(ld,-1);
A = A + diag(ud,1);

% Computing necessary RHS
rhs = A*sol_act;

%% Solving the system using Thomas' algorithm

a = ud;
b = [0 ld];
c = rhs;

for j=2:numel(d)
    d(j) = d(j) - b(j)/d(j-1)*a(j-1);
    c(j) = c(j) - b(j)/d(j-1)*c(j-1);
end 

x = zeros(numel(d),1);
x(numel(d)) = c(numel(d))/d(numel(d));

for j=numel(d)-1:-1:1
    x(j) = (c(j) - a(j)*x(j+1))/d(j);
end