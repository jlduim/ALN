function [S] = espectro(A, tol)

[Q_inicial,R_inicial]=qr_GSM(A); // inicia Q e R como as matrizes obtidas na decomposição de A pelo Método de Gram Schmidt Modificado
A_prox = R_inicial * Q_inicial; // calcula a próxima matriz A, de acordo com o algoritmo implementado
e = norm(diag(A_prox) - diag(A), %inf); // inicia e como a norma infinito

while e > tol do // ao final do laço, a tolerância passada foi atingida
    A = A_prox; // atualiza A para prosseguir obtendo a próxima matriz A
    [Q_prox,R_prox]=qr_GSM(A); // obtém as próximas matrizes Q e R obtidas na decomposição da atual matriz A pelo Método de Gram Schmidt Modificado
    A_prox = R_prox * Q_prox; // obtém a próxima matriz A, de acordo com o algoritmo implementado
    e = norm(diag(A_prox) - diag(A), %inf); // calcula a norma infinito da diferença entre o espectro atual e o obtido na iteração anterior
end

S = diag(A_prox); // armazena o espectro a ser retornado

endfunction
