function [Q]=constroi_Q_House(U)

m = size(U,1); // armazena o número de linhas de U
n = size(U,2); // armazena o número de colunas de U
u_1 = U(:,1); // armazena a primeira coluna de U
Q = eye(m,m) - 2 * u_1 * u_1'; // inicia Q como a matriz ortogonal associada a u_1 obtida pelo método 

for k=2:n // ao final do laço, a matriz Q foi corretamente construída
    u_k = U(:,k); // armazena a k-ésima coluna de U
    Q_k = eye(m,m) - 2 * u_k * u_k'; // calcula a matriz ortogonal associada a u_k obtida pelo método
    Q = Q * Q_k; // atualiza Q multiplicando iterativamente à direita as novas matrizes Q_k obtidas
end

endfunction
