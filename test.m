% Generate a random matrix of size nxn
n=6;
m=9;
A=randi([-n,m],m,n);
% C1=rand(n)
A(A==0)=n;
rank(A);
b=randi([-m,m],m,1);

[data,info] = OneNormLP3033954135(A,b);
[y,fval] = linprog(b, [], [], A', zeros(size(A,2),1), -ones(size(A,1),1), ones(size(A,1),1));
disp(y);
disp(fval);


%degenerate case
A = [1,0;0,1;0,0];
b = [1;1;1];
[data,info] = OneNormLP3033954135(A,b);


