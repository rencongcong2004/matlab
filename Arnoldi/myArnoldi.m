function[Q,H,normv] = myArnoldi(A,n,v)



if iscolumn(A)
    m = length(A);
    A=diag(A);
else
    m = size(A, 1);
end

if nargin<3
    v=ones(m,1);
end

normv=norm(v);

Q = zeros(m,n);
H = zeros(n+1,n);
v = v./norm(v); 
Q(:,1) = v;
for k = 1:n
    w = A*Q(:,k);
    for i = 1:k
        v = Q(:,i);
        H(i,k) = w'*v;
        w = w - H(i,k)*v;
    end
    H(k+1,k) = norm(w,2);
    Q(:,k+1) = w / H(k+1,k);
    if H(k+1,k) < 1e-12
        Q = Q(:,1:k); 
        H = H(1:k,1:k);
        break
    end
end