function [xk1,e,k,rk]=Metodo_Gauss_Seidel_1(A,b,x0,E,M,normtipo)

D = diag(diag(A));
L = tril(A,-1)
U = triu(A,1)
invLMD = inv(L + D);
MG = (-1) * invLMD * U;
cG = invLMD * b;
xk = x0;
xk1 = MG * xk + cG;
e = norm(xk1 - xk, normtipo);
k = 1;

while e >= E,
    xk = xk1;
    xk1 = MG * xk + cG;
    e = norm(xk1 - xk, normtipo);
    k = k + 1;
    if k > M then
        disp("Número máximo de iterações atingido.");
        abort
    end
end

rk = norm(b - A * xk1, normtipo);
    
endfunction
