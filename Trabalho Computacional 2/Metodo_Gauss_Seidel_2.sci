function [x]=Resolve_Sistema_Lx(L,b)

[n]=size(L,1);
x=zeros(n,1);
x(1)=b(1)/L(1,1);
for i=2:n
    x(i)=(b(i)-L(i,1:i)*x(1:i))/L(i,i);
end

endfunction


function [xk1,e,k,rk]=Metodo_Gauss_Seidel_2(A,b,x0,E,M,normtipo)

D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
xk = x0;
xk1 = Resolve_Sistema_Lx(L+D,b-U*xk);
e = norm(xk1 - xk, normtipo);
k = 1;

while e >= E,
    xk = xk1;
    xk1 = Resolve_Sistema_Lx(L+D,b-U*xk);
    e = norm(xk1 - xk, normtipo);
    k = k + 1;
    if k > M then
        disp("Número máximo de iterações atingido.");
        abort
    end
end

rk = norm(b - A * xk1, normtipo);
    
endfunction
