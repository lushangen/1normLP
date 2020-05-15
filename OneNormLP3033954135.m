function[data,info] = OneNormLP3033954135(A,b)

info = infoClass;
data = dataClass;

[m,n] = size(A);
flag = false;
B = [1:n];
B_bar = setdiff([1:m],B);

if ~check_degeneracy(A,b)
    disp('This LP is degenerate.');
    info.run = 'Failure';
    info.msg = 'This LP is degenerate.';
    return
end

x = (A(B,:))\b(B);

if ~check_degeneracy_two(A,b,x)
    disp('This LP is degenerate.');
    info.run = 'Failure';
    info.msg = 'This LP is degenerate.';
    return
end

initx = x;
h = A*x - b; 
h(B_bar) = A(B_bar,:)*x - b(B_bar);
%preallocate y
y = ones(m,1);
y(B_bar) = sign(h(B_bar));
y(B) = -inv(A(B,:))' * A(B_bar,:)'* y(B_bar);

while ~flag
    for i = 1:n
        temp2 = y(B);
        if abs(temp2(i)) > 1 
            j = B(i);
            e = zeros(1,n);
            e(i) = 1;
            t_B_bar = -sign(y(j))*y(B_bar).*(A(B_bar,:)*(A(B,:)\e'));
            temp = abs(h(B_bar))./t_B_bar;
            temp(t_B_bar <= 0) = Inf;
            [~,r] = min(temp);
            B = sort([setdiff(B,j),B_bar(r)]);
            B_bar = setdiff([1:m],B);
            flag = false;
            break;
        end
        flag = true;
    end
    x = (A(B,:))\b(B);
    
    if ~check_degeneracy_two(A,b,x)
        disp('This LP is degenerate.');
        info.run = 'Failure';
        info.msg = 'This LP is degenerate.';
        return
    end
    
    h = A*x - b; 
    h(B_bar) = A(B_bar,:)*x - b(B_bar);
    y(B_bar) = sign(h(B_bar));
    y(B) = -inv(A(B,:))' * A(B_bar,:)'* y(B_bar);
end

obj_val = norm(A(B_bar,:)*x - b(B_bar), 1);


info.run = 'Success';
info.msg = 'No issues';
data.Opt = obj_val;
data.Initx = initx;
data.Optx = x;
data.Opty = y(B);

disp(info.run);
disp(info.msg);
disp('The optimal was');
disp(data.Opt);
disp('With initial x');
disp(data.Initx);
disp('With optimal x');
disp(data.Optx);
disp('With optimal y');
disp(data.Opty);

end

function result = check_degeneracy(A,b)
result = true;

[m,n] = size(A);

B = combnk(1:m,n);
for i = 1:size(B,1)
    if rank(A(B(i,:),:)) ~= n
        result = false;
        return;
    end
end
end

function result = check_degeneracy_two(A,b,x)
result = true;

[m,n] = size(A);

for i = n+1:m
    B = combnk(1:m,i);
    for j = 1:size(B,1)
        if A(B(j,:),:)*x == b(B(j,:))
            result = false;
            return
        end
    end
end

end