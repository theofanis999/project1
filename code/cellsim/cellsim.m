%% introductory cleanup
clear
clc
%% initialization

N = 50;
T = 100;
A = zeros(N,T);
B = zeros(N,T);
db = 0.05;
%% temporal loop
for t=1:T
    %% nodal loop
    for i = 1:N
        r = rand(1);
        A(i,t) =r;
        B(i,t) = B(i,t) + t*db*A(i,t);
    end
end

%% results presentation

disp(B);
disp(A);
