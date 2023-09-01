# Avaliação de estimadores de máxima verossimilhança para os parâmetros de escala e forma da distribuição de Pareto generalizada

Neste trabalho, busca-se obter estimativas de MV dos parâmetros de escala e forma, $\sigma$ e $\xi$, respectivamente, da GPD, por meio de simulação de Monte Carlo. Buscamos ainda corrigir estas estimativas utilizando a amostragem pelo método *bootstrap*.

Os estimadores de MV são viesados em relação ao verdadeiro valor do parâmetros, sendo portanto influenciados pelo tamanho da amostra. Diferentes estratégias podem ser utilizadas para verificar as caracteristicas de um estimador. Simulações de Monte Carlo são úteis pois estas geram uma distribuição empírica do estimador, permitindo a avaliação do viés, variância e erro quadrático médio deste em relação a um valor conhecido do parâmetro. Outras técnicas, como a amostragem por *bootstrap*, podem ser empregadas em conjunto com o processo de simulação de Monte Carlo para a diminuição do viés do estimador.

A GPD é definida por pelos parâmetros de escala ($\sigma$) e de forma ($\xi$), com $\sigma>0$ e $-\infty <\xi< \infty$. Se $X \sim GPD(\sigma, \xi)$, sua função de distribuição acumulada é:

```math
\begin{equation*}
  F(x;\sigma, \xi) =
    \begin{cases}
      1 - \left( 1 - \frac{\xi x}{\sigma} \right)^{\frac{1}{\xi}} & \xi \ne 0\\
      1 - e^{\frac{-x}{\sigma}} & \xi=0
    \end{cases}       
\end{equation*}
```
com suporte $x>0$ para $\xi \le 0$ e $0 < x < \frac{\sigma}{\xi}$ para $\xi>0$. Além disso, 

```math
\begin{equation*}
  E(X) = \frac{\sigma}{1-\xi} \text{, para $\xi < 1$}
\end{equation*}
```

```math
\begin{equation*}
  Var(X) = \frac{\sigma^2}{(1-\xi)^2(1-2\xi)} \text{, para $\xi < 1/2$}
\end{equation*}
```
