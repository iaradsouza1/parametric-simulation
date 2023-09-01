# Evaluation of maximum likelihood estimators for the generalized Pareto distribution (GPD) parameters

In this work, we seek to obtain MLE from the scale and shape parameters, $\sigma$ and $\xi$, respectively, from GPD, by means of Monte Carlo simulation. Also, we seek to correct these estimates using the *bootstrap* method.

The MLE are compared to the parameters true value, being therefore influenced by the size of the sample. Different strategies can be used to verify the characteristics of an estimator. Monte Carlo simulations are useful because they generate an empirical distribution of the estimator, allowing the evaluation of bias, variance, and mean square error in comparison to the parameter value. Other techniques, such as bootstrap sampling, can be used with the Monte Carlo simulation process to decrease the estimator's range.

The GPD is defined by the scale ($\sigma$) and shape ($\xi$) parameters, with $\sigma>0$ and $-\infty <\xi< \infty$. If $X \sim GPD(\sigma, \xi)$, its cumulative distribution function is:

```math
\begin{equation*}
  F(x;\sigma, \xi) =
    \begin{cases}
      1 - \left( 1 - \frac{\xi x}{\sigma} \right)^{\frac{1}{\xi}} & \xi \ne 0\\
      1 - e^{\frac{-x}{\sigma}} & \xi=0
    \end{cases}       
\end{equation*}
```
with support of $x>0$ for $\xi \le 0$ and $0 < x < \frac{\sigma}{\xi}$ for $\xi>0$. Additionally, 

```math
\begin{equation*}
  E(X) = \frac{\sigma}{1-\xi} \text{, for $\xi < 1$}
\end{equation*}
```

```math
\begin{equation*}
  Var(X) = \frac{\sigma^2}{(1-\xi)^2(1-2\xi)} \text{, for $\xi < 1/2$}
\end{equation*}
```
