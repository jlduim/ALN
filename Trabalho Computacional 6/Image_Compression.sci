function [C]=Image_Compression(A,p)

A = double(A); // converte as entradas da matriz A para o tipo double
[U,S,V] = svd(A); // realiza a decomposição SVD da matriz A
r = length(diag(S)); // realiza a contagem de valores singulares da matriz S e armazena na variável r
s = max(1, floor(p * r)); // calcula quantos valores singulares serão mantidos na imagem comprimida, dada a porcentagem p passada
C = U(:,1:s) * S(1:s,1:s) * V'(1:s,:); // obtêm a matriz C correspondente à imagem comprimida
C = iconvert(C,11); // converte as entradas da matriz C para o tipo inteiro com 1 byte (0 a 255)

endfunction
