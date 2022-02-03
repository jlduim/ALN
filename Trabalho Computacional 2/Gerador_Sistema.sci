function [A,b]=Gerador_Sistema(n)

M = rand(n,n);
for i=1:n
    M(i,i) = sum(abs(M(i,1:n))) + 1;
end

A = M;
b = rand(n,1);
    
endfunction
