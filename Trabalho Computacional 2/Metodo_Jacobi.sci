function [xk1,e,k,rk]=Metodo_Jacobi(A,b,x0,E,M,normtipo)

D = diag(diag(A));
LMU = A - D;
invD = diag(1./diag(D));
MJ = (-1) * invD * LMU;
cJ = invD * b;
xk = x0;
xk1 = MJ * xk + cJ;
e = norm(xk1 - xk, normtipo);
k = 1;

while e >= E,
    xk = xk1;
    xk1 = MJ * xk + cJ;
    e = norm(xk1 - xk, normtipo);
    k = k + 1;
    if k > M then
        disp("Número máximo de iterações atingido.");
        abort
    end
end

rk = norm(b - A * xk1, normtipo);
    
endfunction
