function [Q,R]=qr_GSM(A)

m = size(A,1); // armazena o número de linhas de A
n = size(A,2); // armazena o número de colunas de A
Q = eye(m,n); // inicia Q como uma matriz "identidade" de mesmas dimensões de A
R = eye(n,n); // inicia R como uma matriz identidade n x n

for j=1:n // ao final do laço, Q e R foram corretamente determinadas
    a_j = A(:,j); // armazena a j-ésima coluna de A
    v = a_j; // inicia v como a j-ésima coluna de A
    for i=1:j-1 // ao final do laço, v foi ortogonalizado em relação ao subespaço gerado pelas (j-1) primeiras colunas de A
        q_i = Q(:,i); // armazena a i-ésima coluna de Q (já calculada em iterações anteriores)
        r_ij = q_i' * v; // calcula a entrada ij de R
        R(i,j) = r_ij; // insere a entrada ij de R na matriz R
        v = v - r_ij * q_i; // atualiza v como a sua componente ortogonal em relação ao subespaço gerado pela i-ésima coluna de A
    end
    r_jj = norm(v); // armazena a norma de v
    R(j,j) = r_jj; // insere a norma de v na j-ésima entrada da diagonal principal de R
    q_j = v / r_jj; // normaliza v
    Q(:,j) = q_j; // insere o vetor v normalizado na j-ésima coluna de Q
    
end

endfunction
