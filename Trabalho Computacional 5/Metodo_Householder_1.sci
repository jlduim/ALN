function [U,R]=qr_House(A)

m = size(A,1); // armazena o número de linhas de A
n = size(A,2); // armazena o número de colunas de A
U = zeros(m,n); // inicia U como uma matriz nula de mesmas dimensões de A
R = eye(m,n); // inicia R como uma matriz "identidade" de mesmas dimensões de A

for k=1:n // ao final do laço, U e R foram corretamente determinadas
    x = A(k:m,k); // inicia x como os elementos da k-ésima coluna de A excluindo os elementos acima da diagonal principal
    if x(1) < 0 then // verifica se a primeira entrada de x é negativa
        x(1) = x(1) - norm(x); // subtrai a norma de x de sua primeira entrada
    else // verifica se a primeira entrada de x é positiva ou nula
        x(1) = x(1) + norm(x); // soma a norma de x à sua primeira entrada
    end
    u = x / norm(x); // normaliza x
    U(k:m,k) = u; // insere o vetor x normalizado na k-ésima coluna de U excluindo os elementos acima da diagonal principal
    A(k:m,k:n) = A(k:m,k:n) - 2 * u * (u' * A(k:m,k:n)); // atualiza, segundo o método implementado, a submatriz de A obtida excluindo as (k-1) primeiras linhas e colunas de A
end

R = triu(A); // substitui os elementos de A abaixo da diagonal principal por zeros

endfunction
