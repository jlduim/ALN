function [lambda,x1,k,n_erro]=Metodo_Potencia_2(A,x0,epsilon,M)

k = 1; // inicia a contagem de iterações
x0 = x0 / norm(x0, 2); // normalização de x0
x1 = A * x0; // execução da primeira iteração

while k <= M,
    lambda = x0' * x1; // calcula o quociente de Rayleigh
    if lambda < 0 then // verifica se lambda é negativo
        x1 = - x1; // inverte a orientação de x1, mantendo-o com o mesmo sentido de x0
    end
    x1 = x1 / norm(x1, 2); // normalização de x1
    n_erro = norm(x1 - x0, 2); // quantifica o erro entre as 2 últimas iterações
    if n_erro < epsilon then // verifica se a precisão desejada foi atingida
        disp("Convergência atingida.");
        break // interrompe as iterações
    end
    x0 = x1; // atualiza o valor de x0 com a última iteração
    x1 = A * x0; // execução da próxima iteração
    k = k + 1; // incrementa a contagem de iterações em 1 unidade
end

if k > M then // verifica se o número de iterações ultrapassou o limite dado
    disp("Número máximo de iterações atingido.");
end

endfunction
