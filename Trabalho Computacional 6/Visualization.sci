function []=Visualization(A,p)

C = Image_Compression(A,p); // comprime a imagem e armazena a matriz correspondente na variável C
subplot(1,2,1); imshow(A); // plota a imagem original no 1° espaço (à esquerda) reservado pelo subplot
subplot(1,2,2); imshow(C); // plota a imagem comprimida no 2° espaço (à direita) reservado pelo subplot

endfunction
