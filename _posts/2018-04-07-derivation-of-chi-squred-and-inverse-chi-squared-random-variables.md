---
layout: post
title: Derivation of Chi-squred and Inverse Chi-squared Distributions
date: 2018-04-07 04:23
author: baruuum
comments: true
categories: [Chi-squared, Distributions, Inverse Chi-squared, Quant Stuff]
---
In a recent project, I had to simulate multivariate t-distributed random variates. While there are functions which directly sample from multivarite t-distributions, it is convenient (and often faster) to derive these variates from other distributions which are easier to sample from. For example, to sample from a univariate t distribution with $ k $ degrees of freedom, we might first sample $ Z $ a from standard Normal distribution, then sample independently $ S $ from a Chi-squared distribution with $ k $ degrees of freedom, and calculate $ T = Z/\sqrt{S/k} $ which will be a sample from the desired distribution. Anyways, this brought me to inverse Chi-squared variables, which I had never encountered before. It turned out that the derivation of the inverse Chi-squared distribution is quite straightforward. So, I use the spare time while my models are running to post about them.

First, we start from a random variable $ Z $ that follows a standard Normal distribution. It is well-known that $ X=Z^2 $ follows a chi-squared distribution. Thus, let us derive the distribution of $ X $. Denote by $ \Phi $ the standard Normal CDF. For $ x > 0 $, $ X $ has CDF

$$ \begin{aligned}
F_X(x) &=\Pr[X \le x] = \Pr[Z^2 \le x] = \Pr[\vert Z\vert  \le \sqrt{x}] = \Pr[ -\sqrt x \le Z \le \sqrt x] \\
&= \int_{-\sqrt x}^{\sqrt x} \frac{1}{\sqrt{2\pi}} e^{-\frac{z^2}{2}}dz = 2 \int_{0}^{\sqrt x} \frac{1}{\sqrt{2\pi}} e^{-\frac{z^2}{2}}dz.
\end{aligned}$$

So, it follows that the pdf of $ X $ is

$$ \begin{aligned}
f(x) = \frac{d}{dx} F_X(x) = \frac{d}{dx}[\Phi(\sqrt{x}) - \Phi(0)] = 2\phi(\sqrt{x}) \frac{d\sqrt x}{dx} = \frac{1}{2^{1/2}\Gamma(1/2)} x^{-1/2} e^{-x/2}, \quad x>0
\end{aligned}$$

where $ \Gamma(\cdot) $ is gamma function and $ \Gamma(1/2)=\sqrt \pi $. Note that this distribution is a special case of the Gamma distribution with shape and rate parameter equal to 1/2 and 1/2, respectively. In general, the Gamma distribution with shape parameter $ \alpha $ and rate parameter $ \beta $ is given as

$$  g(w; \alpha,\beta) = \frac{\beta^\alpha}{\Gamma(\alpha)} w^{\alpha - 1} e^{-\beta w}, \quad w>0. $$

Plugging in $ \alpha = \beta = 1/2 $ gives thus the Chi-squared distribution. It will be convenient to derive the moment generating function (MGF) of the Gamma distribution:

$$  M_W(t) = E[e^{tw}] = \int_{0}^\infty \frac{\beta^\alpha}{\Gamma(\alpha)} w^{\alpha - 1} e^{-(\beta -t)w} dx = \left(1-\frac{t}{\beta}\right)^{-\alpha} $$

from this it follows that the MGF of a Chi-squared variable is $ M_X(t) = (1-2t)^{-1/2} $. Now, consider the sum of $ \nu $ independently distributed Chi-squared random variables. By independence, the MGF of $ Y = \sum_{j=1}^\nu X_j $ is

$$  M_Y(t) = \prod_{j=1}^\nu M_{X_j}(t) = \prod_{j=1}^\nu M_X(t) = (1-2t)^{-\nu/2} $$

which is the MGF of a Gamma distribution with shape parameter $ \nu/2 $ and rate parameter $ 1/2 $. Hence, the distribution of $ Y $ is given as

$$  f_Y(y ; \nu) = \frac{1}{2^{\nu/2}\Gamma(\nu/2)} y^{\nu/2 - 1} e^{-y/2}, \quad y>0, $$

which is known as the pdf of the Chi-squared distribution with $ \nu $ degrees of freedom. It is immediately clear that what we've so far called *the* Chi-squared distribution is a Chi-squared distribution with one degree of freedom.

Finally, we derive the distribution of $ V = 1/Y $, where $ Y \sim \text{Chi-sqaured}(\nu) $. As $ Y $ has support for strictly positive values, the following holds

$$  F_V(v) = \Pr[V \le v ] = \Pr[1/Y \le v] = \Pr[Y > 1/v] = 1- F_Y(1/v) $$

Thus, the pdf of $ V $ is given by the derivative

$$\begin{aligned}
f_V(v ; \nu) &= \frac{d}{dv}\Big[1-F_Y(1/v)\Big] = -\frac{d}{dv}F_Y(1/v) = \frac{1}{v^2} f_Y(1/y) \\
&=\frac{(1/2)^{\nu/2}}{\Gamma(\nu/2)} v^{-\nu/2 - 1} e^{-\frac{1}{2v}}.
\end{aligned}.$$

The distribution with pdf $ f_V(\cdot\, ;\, \nu) $ is known as the inverse Chi-squared distribution with $ \nu $ degrees of freedom. 
