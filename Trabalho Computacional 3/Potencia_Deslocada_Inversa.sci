function [lambda,x1,k,n_erro]=Potencia_Deslocada_Inversa(A,x0,epsilon,alfa,M)

k = 1; // inicia a contagem de iterações
x0 = x0 / norm(x0, 2); // normalização de x0
n = length(x0); // armazena a dimensão de x0 para criar uma matriz identidade com as dimensões corretas

while k <= M,
    x1 = Gaussian_Elimination_4(A - alfa * eye(n, n), x0); // resolve o sistema (A - alfa * I) * x1 = x0 na variável x1 a fim de iterar inversamente
    x1 = x1 / norm(x1, 2); // normalização de x1
    lambda = x1' * A * x1; // calcula o quociente de Rayleigh
    if lambda - alfa < 0 then // verifica se (lambda - alfa) é negativo // alteração feita pois os autovalores de A foram deslocados de alfa ao se fazer A - alfa * I
        x1 = - x1; // inverte a orientação de x1, mantendo-o com o mesmo sentido de x0
    end
    n_erro = norm(x1 - x0, 2); // quantifica o erro entre as 2 últimas iterações
    if n_erro < epsilon then // verifica se a precisão desejada foi atingida
        disp("Convergência atingida.");
        break // interrompe as iterações
    end
    x0 = x1; // atualiza o valor de x0 com a última iteração
    k = k + 1; // incrementa a contagem de iterações em 1 unidade
end

if k > M then // verifica se o número de iterações ultrapassou o limite dado
    disp("Número máximo de iterações atingido.");
end

endfunction
