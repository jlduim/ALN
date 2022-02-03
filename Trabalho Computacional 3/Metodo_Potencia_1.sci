function [lambda,x1,k,n_erro]=Metodo_Potencia_1(A,x0,epsilon,M)

k = 1; // inicia a contagem de iterações
[m_x0, i_x0] = max(abs(x0)); // max(abs(x0)) devolve o módulo da coordenada de maior módulo (m_x0) e o index dessa coordenada (i_x0)
x0 = x0 / x0(i_x0); // x0(i_x0) é a coordenada de maior módulo de x0, utilizada para normalizar o mesmo
x1 = A * x0; // execução da primeira iteração

while k <= M,
    [m_x1, i_x1] = max(abs(x1)); // max(abs(x1)) devolve o módulo da coordenada de maior módulo (m_x1) e o index dessa coordenada (i_x1)
    lambda = x1(i_x1); // x1(i_x1) é a coordenada de maior módulo de x1
    x1 = x1 / lambda; // normalização de x1
    n_erro = norm(x1 - x0, %inf); // quantifica o erro entre as 2 últimas iterações
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
