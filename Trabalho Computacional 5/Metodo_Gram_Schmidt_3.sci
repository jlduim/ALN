function [Q,R,P]=qr_GSP(A)

m = size(A,1); // armazena o número de linhas de A
n = size(A,2); // armazena o número de colunas de A
Q = eye(m,n); // inicia Q como uma matriz "identidade" de mesmas dimensões de A
R = eye(n,n); // inicia R como uma matriz identidade n x n
P = eye(n,n); // inicia P como uma matriz identidade n x n
B = A; // inicia B como sendo igual à matriz A

aux = []; // inicia um vetor auxiliar
for k=1:n // ao final do laço, as normas de todas as colunas de A foram armazenadas
    aux(k) = norm(A(:,k)); // armazena a norma da k-ésima coluna de A na k-ésima entrada do vetor auxiliar
end
[p,q] = max(aux); // obtém, respectivamente, a maior dentre as normas das colunas de A e o número da coluna de maior norma de A
if q ~= 1 then // verifica se a coluna de maior norma de A não é a primeira coluna
    A(:,[1 q]) = A(:,[q 1]); // faz a troca entre a primeira coluna de A e a coluna de maior norma de A
    P(:,[1 q]) = P(:,[q 1]); // faz as mesmas trocas de colunas efetuadas em A de forma que a fatoração obtida ao final corresponda à matriz A original multiplicada por P à direita
end
a_1 = A(:,1); // armazena a primeira coluna de A
R(1,1) = norm(a_1); // insere a norma da primeira coluna de A na j-ésima entrada da diagonal principal de R
Q(:,1) = a_1 / norm(a_1); // insere a primeira coluna de A normalizada na primeira coluna de Q

for j=2:n // ao final do laço, P, Q e R foram corretamente determinadas
    for t=j:n // ao final do laço, as colunas de A ainda não processadas foram ortogonalizadas em relação ao subespaço gerado pela coluna de A sendo sendo processada
        A(:,t) = A(:,t) - ((Q(:,j-1)' * A(:,t)) / (Q(:,j-1)' * Q(:,j-1))) * Q(:,j-1); // ortogonaliza a t-ésima coluna de A em relação ao subespaço gerado pela coluna de A sendo processada
    end
    aux = []; // limpa os dados do vetor auxiliar
    for k=j:n // ao final do laço, as normas de todas as colunas de A ainda não processadas foram armazenadas
        aux(k-j+1) = norm(A(:,k)); // armazena a norma da k-ésima coluna de A na (k-j+1)-ésima entrada do vetor auxiliar
    end
    [p,q] = max(aux); // obtém, respectivamente, a maior dentre as normas das colunas de A e o número da coluna de maior norma de A
    if q ~= 1 then // verifica se a coluna de maior norma de A não é a j-ésima coluna
        A(:,[j q+j-1]) = A(:,[q+j-1 j]); // faz a troca entre a j-ésima coluna de A e a coluna de maior norma de A
        P(:,[j q+j-1]) = P(:,[q+j-1 j]); // faz as mesmas trocas de colunas efetuadas em A de forma que a fatoração obtida ao final corresponda à matriz A original multiplicada por P à direita
    end
    a_j = A(:,j); // armazena a j-ésima coluna de A
    v = a_j; // inicia v como a j-ésima coluna de A
    for i=1:j-1 // ao final do laço, v foi ortogonalizado em relação ao subespaço gerado pelas (j-1) primeiras colunas de Q
        q_i = Q(:,i); // armazena a i-ésima coluna de Q
        prod_esc = q_i' * v; // obtém o produto escalar entre a i-ésima coluna de Q e o vetor v
        aux = []; // limpa os dados do vetor auxiliar
        for z=1:n // ao final do laço, o produto entre a i-ésima linha de Q e cada coluna de B foram armazenados
            aux(z) = Q(:,i)' * B(:,z); // armazena o produto entre a i-ésima linha de Q e a z-ésima coluna de B na z-ésima entrada do vetor auxiliar
        end
        R(i,j) = aux' * P(:,j); // insere o produto escalar entre o vator auxiliar e a j-ésima coluna de P na entrada ij de R
        v = v - prod_esc * q_i; // ortogonaliza v em relação ao subespaço gerado pela i-ésima coluna de Q
    end
    r_jj = norm(v); // armazena a norma de v
    R(j,j) = r_jj; // insere a norma de v na j-ésima entrada da diagonal principal de R
    q_j = v / r_jj; // normaliza v
    Q(:,j) = q_j; // insere o vetor v normalizado na j-ésima coluna de Q
end

endfunction
