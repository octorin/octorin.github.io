---
layout: post
title: Gibbs Sampling
date: 2018-03-15 06:17
author: baruuum
comments: true
categories: [Gibbs Sampling, MCMC, Monte Carlo Methods, Quant Stuff]
---
In this post, I dig a little bit into the theory behind the Gibbs sampler. We will see that by repeatedly sampling from conditional densities we can sample from the joint density. When I first came across this result, it seemed too simple to be true. Yet, as simple as the result seems, the derivation of it is not trivial. Here, again, we will not go too deep into the statistical theory but only try to get a feeling of why the algorithm works.

Let $ \theta  $ be the parameter vector (here treated as a random variable) and suppose that it can be partitioned into $ k  $ subsets, so that $ \theta= (\theta_1,\theta_2,...,\theta_k)  $. Note that each of the $ \theta_j  $'s might be multidimensional as well. Denote the support of $ \theta_j  $ by $ \Theta_j  $ and let $ \Theta  $ be the support of $ \theta  $. Lastly, let $ p  $ be the joint density of $ \theta  $, $ p_j(\theta_j\vert\cdot)  $ the <em>conditional</em> density of $ \theta_j  $, and $ p_j(\cdot)  $ be the marginal density of $ \theta_j  $. With these definitions, the Gibbs sampler can be defined as follows:

- For $ t=1,2,...,T  $,
1. sample $ \theta_1^{t+1} \sim p_1(\theta_1\vert\theta_2^t,\theta_3^t,...,\theta_k^t)  $
2. sample $ \theta_2^{t+1} \sim p_2(\theta_2\vert\theta_1^{t+1},\theta_3^t,...,\theta_k^t)  $
3. ...
4. sample $ \theta_k^{t+1} \sim p_k(\theta_k\vert\theta_1^{t+1},\theta_2^{t+1},...,\theta_{k-1}^{t+1})  $
5. set $ \theta^{(t+1)} = (\theta_1^{t+1}, \theta_2^{t+1},...,\theta_k^{t+1})^\top  $


It turns out that this simple algorithm generates a Markov chain with unique stationary distribution $ p(\theta)  $. MAGIC.

<h3>Joint and Conditional Distributions</h3>

To illustrate why the Gibbs sampler works, we first focus on the case where $ \theta=(\theta_1,\theta_2)  $. We establish the result that the conditional densities completely characterize the joint density, given that the integral $ \int [p_2(x\vert\theta_1)/p_1(\theta_1\vert x)]dx  $ exists. With this assumption we have

$$  \int\frac{p_2(x\vert\theta_1)}{p_1(\theta_1\vert x)}dx = \int \frac{p(x,\theta_1)/p_1(\theta_1)}{p(x,\theta_1)/p_2(x)}dx = \frac{1}{p_1(\theta_1)}\int p_2(x)dx = \frac{1}{p_1(\theta_1)}.  $$

Therefore,

$$  p(\theta_1,\theta_2) = p_2(\theta_2\vert\theta_1)p_1(\theta_1) = \frac{p_2(\theta_2\vert\theta_1)}{\int [p_2(x\vert\theta_1)/p_1(\theta_1\vert x)]dx}.  $$

Thus, the joint density $ p(\theta_1,\theta_2)  $ can be expressed as a function of the conditional densities $ p_1(\theta_1\vert\theta_2)  $ and $ p_2(\theta_2\vert\theta_1)  $.

Note that the following condition is sufficient to ensure that the integral exists:

$$  p_j(\theta_j)>0 \text{ for all }j=1,2,..,k \implies p(\theta)>0.  $$

This condition, called the <em>positivity condition</em>, states that if the events $ \theta_1  $ and $ \theta_2  $ can occur separately, then they can also occur together. In the bivariate example from above, it implies that $ p_1(\theta_1\vert\theta_2), p_2(\theta_2\vert\theta_1)>0  $, whenever $ p_1(\theta_1), p_2(\theta_2)>0  $ as $ p(\theta_1,\theta_2) = p_1(\theta_1\vert\theta_2)p_2(\theta_2)= p_2(\theta_2\vert\theta_1) p(\theta_1)  $.

The generalization of this result will lead us to the famous <em>Hammersley-Clifford Theorem</em>. Suppose we want to sample from the joint density $ p(\theta)  $. Let $ \omega, \xi \sim p  $. Assume that if $ p_i(\omega_i)>0, p_j(\omega_j)>0  $ for all $ i,j  $, then $ p(\omega,\xi)>0  $ (the positivity condition). First, consider

$$  p(\omega) = p(\omega_n\vert\omega_2,...,\omega_{n-1})p(\omega_1,...,\omega_{n-1})  $$

Next, we note that $ p(\xi_n\vert\omega_1,...,\omega_{n-1})p(\omega_1,...,\omega_{n-1}) = p(\omega_1,...,\omega_{n-1}, \xi_n)  $, so that

$$  p(\omega) = \frac{p(\omega_n\vert\omega_1,...,\omega_{n-1})}{p(\xi_n\vert\omega_1,...,\omega_{n-1})}p(\omega_1,...,\omega_{n-1},\xi_n).  $$

Using the same trick on $ p(\omega_1,...,\omega_{n-1},\xi_n)  $ with respect to $ \omega_{n-1}  $ and introducing $ \xi_{n-1}  $ leads to

$$ \begin{aligned}
p(\omega) &= \frac{p(\omega_n\vert\omega_1,...,\omega_{n-1})}{p(\xi_n\vert\omega_1,...,\omega_{n-1})}p(\omega_1,...,\omega_{n-1},\xi_n) \\
&= \frac{p(\omega_n\vert\omega_1,...,\omega_{n-1})}{p(\xi_n\vert\omega_1,...,\omega_{n-1})}p(\omega_{n-1}\vert\omega_1,...,\omega_{n-2},\xi_n)p(\omega_1,...,\omega_{n-2},\xi_{n}) \\
&=\frac{p(\omega_n\vert\omega_1,...,\omega_{n-1})}{p(\xi_n\vert\omega_1,...,\omega_{n-1})}\frac{p(\omega_{n-1}\vert\omega_1,...,\omega_{n-2},\xi_n)}{p(\xi_{n-1}\vert\omega_1,...,\omega_{n-2},\xi_n)}p(\omega_1,...,\omega_{n-2},\xi_{n-1},\xi_n) \\
&\quad \vdots \\
&=\prod_{i=1}^n\frac{p(\omega_i\vert\omega_1,...,\omega_{i-1},\xi_{i+1},...,\xi_n)}{p(\xi_i\vert\omega_1,...,\omega_{i-1},\xi_{i+1},...,\xi_n)}p(\xi_1,...,\xi_n).
\end{aligned}
$$

Thus, we have

$$  \frac{p(\omega)}{p(\xi)} = \prod_{i=1}^n\frac{p(\omega_i\vert\omega_1,...,\omega_{i-1},\xi_{i+1},...,\xi_n)}{p(\xi_i\vert\omega_1,...,\omega_{i-1},\xi_{i+1},...,\xi_n)},  $$

or

$$  p(\omega) \propto \prod_{i=1}^n\frac{p(\omega_i\vert\omega_1,...,\omega_{i-1},\xi_{i+1},...,\xi_n)}{p(\xi_i\vert\omega_1,...,\omega_{i-1},\xi_{i+1},...,\xi_n)},  $$

which is a charactarization of the joint density by conditionals.

A concrete example might help. Suppose we want to sample from a bivariate joint distribution of $ \theta_1  $ and $ \theta_2  $, knowing that the conditional distributions are

$$  \theta_i\vert\theta_j \sim \text{Normal}(\rho\theta_j, (1-\rho^2)).  $$

The Hammersley-Clifford Theorem tells us that

$$ \begin{aligned}
p(\theta_1,\theta_2) &\propto \frac{p(\theta_1\vert x_2)p(\theta_2\vert\theta_1)}{p(x_1\vert x_2)p(x_2\vert\theta_1)} \\
&\propto \frac{\exp\left[-\frac{1}{2(1-\rho^2)}(\theta_1 - \rho x_2)^2\right] \exp\left[-\frac{1}{2(1-\rho^2)}(\theta_2 - \rho \theta_1)^2\right] }{\exp\left[-\frac{1}{2(1-\rho^2)}(x_1 - \rho x_2)^2\right]\exp\left[-\frac{1}{2(1-\rho^2)}(x_2 - \rho \theta_1)^2\right]} \\
&\propto \exp\left[-\frac{\theta_1^2 + \theta_2^2 - 2\rho \theta_1\theta_2}{2(1-\rho^2)} \right],
\end{aligned} $$

which is the kernel of a bivariate Normal distribution with mean vector $ \mu=[0,0]  $ and covariance matrix $$ \Sigma=\begin{bmatrix} 1&\rho\\ \rho& 1\end{bmatrix}.$$

Before we return to the Gibbs sampler, a note of caution is needed. While it is true that the conditional distributions exist if the joint distribution exists, the reverse of the statement is not necessarily true. In other words, there are situations in which two (or more) valid conditional distributions lead to a joint distribution which is undefined. Yet, the Gibbs sampler, will just generate samples from whatever the conditional distributions suggest (it is just an algorithm!). Hence, care is needed when applying the Gibbs sampler in practice.

<h3>Markov Chain Generated by the Gibbs Sampler</h3>

Now, consider the transition kernel of the Gibbs sampler, i.e.,

$$  P(\theta^t, \theta^{t+1}) = p_1(\theta_1^{t+1}\vert\theta_2^t)\cdot p_2(\theta_2^{t+1}\vert\theta_1^{t+1}).  $$

Suppose the sampler is at $ \theta^t  $ and consider the probability that it visits a subset $ A  $ of $ \Theta  $ at iteration $ t+1  $:

$$  \int_A P(\theta^t, \theta^{t+1}) d\theta^{t+1} = \int_A p_1(\theta_1^{t+1}\vert\theta_2^t)\cdot p_2(\theta_2^{t+1}\vert\theta_1^{t+1})d\theta^{t+1}.  $$

If the positivity condition holds, then this transition probability will be positive for any $ A\subset \Theta  $, so that, at each step, the Gibbs sampler can visit any subset of $ A  $ of $ \Theta  $. Thus, under positivity, the Markov chain generated by the Gibbs sampler is irreducible (all states communicate), strongly reccurent (any state has a non-zero probability of being visited in one step from any other state), and aperiodic (there is a non-zero probability to stay at a state, implying that the period of the communication class is one). It follows that it has a unique stationary distribution.

Some subtleties have to be added if $ \theta  $ is a continuous random variable as the probability of moving to any point (i.e., a singleton set) in $ \Theta  $ will be always zero. So, formulating the transition probabilities in terms of "arbitrary" subsets of sets $ \Theta  $ is, strictly speaking, wrong. To be rigorous, we should have said that the Markov Chain generated by the Gibbs sampler is able to jump to every set $ A \in \mathcal B(\Theta)  $ with positive probability measure, where $ \mathcal B(\Theta)  $ is the Borel sigma-algebra on $ \Theta  $.

It remains to show that the the stationary distribution to which the Gibbs sampler converges coincides with $ p  $, the joint distribution of $ \theta_1  $ and $ \theta_2  $ from which we want to sample. Suppose that we start the Gibbs sampler at the target distribution $ p  $, so that $ \theta^{t-1} \sim p  $. We want to show that the probability that the Gibbs sampler will visit an arbitrary subset $ A  $ with positive probability in the next step is equal to the probability under the joint distribution $ p  $. This would imply that the Gibbs sampler, once reaching the distribution $ p  $, will stay there. Thus, we have to show that $ \Pr[\theta^t \in A] = \int_A p(\theta^{t})d\theta^{t}  $ for all measurable $ A  $, assuming that $ \theta^{t-1}\sim p  $. So, let us calculate this probability.

$$\begin{aligned}
\Pr[\theta^t \in A] &= \int_A \int_\Theta P(\theta^{t-1},\theta^t) p(\theta^{t-1}) d\theta^{t-1}d\theta^t \\
&=\int_\Theta \int_\Theta \mathbb I(\theta^t \in A) P(\theta^{t-1},\theta^t) p(\theta^{t-1}) d\theta^{t-1}d\theta^t.
\end{aligned} $$

But $ P(\theta^{t-1},\theta^t) = p_1(\theta_1^{t}\vert\theta_2^{t-1}) p_2(\theta_2^{t}\vert\theta_1^{t})  $, so that

$$  \Pr[\theta^t \in A]=\int_\Theta \int_\Theta \mathbb I(\theta^t \in A)p_2(\theta_2^{t}\vert\theta_1^{t}) p_1(\theta_1^{t}\vert\theta_2^{t-1}) p(\theta^{t-1}) d\theta^{t-1}d\theta^t.  $$

Note that $ \theta_1^{t-1}  $ enters the expression only through $ p(\theta^{t-1})  $. Thus, we might integrate it out noting that $ \int_{\Theta_1} p(\theta^{t-1})d\theta_1^{t-1} = p_2(\theta_2^{t-1})  $. Since $ p_1(\theta_1^{t}\vert\theta_2^{t-1})p_2(\theta_2^{t-1})= p(\theta_1^t, \theta_2^{t-1})  $, we have

$$\begin{aligned}
\Pr[\theta^t \in A]&=\int_{\Theta_2}\int_{\Theta_1}\int_{\Theta_2}\int_{\Theta_1} \mathbb I(\theta^t \in A)p_2(\theta_2^{t}\vert\theta_1^{t}) p_1(\theta_1^{t}\vert\theta_2^{t-1}) p(\theta^{t-1}) d\theta_1^{t-1}d\theta_2^{t-1}d\theta_1^{t}d\theta_2^{t} \\
&=\int_{\Theta_2}\int_{\Theta_1}\int_{\Theta_2} \mathbb I(\theta^t \in A)p_2(\theta_2^{t}\vert\theta_1^{t}) p_1(\theta_1^{t}\vert\theta_2^{t-1})\left( \int_{\Theta_1}p(\theta^{t-1}) d\theta_1^{t-1}\right)d\theta_2^{t-1}d\theta_1^{t}d\theta_2^{t} \\
&=\int_{\Theta_2}\int_{\Theta_1}\int_{\Theta_2} \mathbb I(\theta^t \in A)p_2(\theta_2^{t}\vert\theta_1^{t}) p_1(\theta_1^{t}\vert\theta_2^{t-1})p_2(\theta_2^{t-1})d\theta_2^{t-1}d\theta_1^{t}d\theta_2^{t} \\
&=\int_{\Theta_2}\int_{\Theta_1}\int_{\Theta_2} \mathbb I(\theta^t \in A)p_2(\theta_2^{t}\vert\theta_1^{t}) p(\theta_1^t,\theta_2^{t-1})d\theta_2^{t-1}d\theta_1^{t}d\theta_2^{t} \\
&=\int_{\Theta_2}\int_{\Theta_1} \mathbb I(\theta^t \in A)p_2(\theta_2^{t}\vert\theta_1^{t})\left(\int_{\Theta_2} p(\theta_1^t,\theta_2^{t-1})d\theta_2^{t-1}\right)d\theta_1^{t}d\theta_2^{t}\\
&= \int_{\Theta_2}\int_{\Theta_1} \mathbb I(\theta^t \in A)p_2(\theta_2^{t}\vert\theta_1^{t})p_1(\theta^t)d\theta_1^{t}d\theta_2^{t}\\
&= \int_{\Theta_2}\int_{\Theta_1} \mathbb I(\theta^t \in A)p(\theta_1^t,\theta_2^{t})d\theta_1^{t}d\theta_2^{t}\\
&=\int_\Theta \mathbb I(\theta^t \in A) p(\theta)^t d\theta^t \\
&= \int_A p(\theta^t)d\theta^t
\end{aligned} $$

as desired. Hence, the Markov chain generated by the Gibbs sampler has the joint density as the unique stationary distribution.

The Gibbs sampler is very useful in situations where sampling form the joint density is intractable, but sampling from the conditionals is easy. A situation like this arises naturally when estimating mixture models, which will be our example below. Suppose that we want to model a distribution as the mixture or Normally distributed random variables, where the number of latent classes and the class-specific variances are known in advance. We are interested in the posterior distribution of $ (\lambda, \mu)  $, where $ \lambda=(\lambda_1,\lambda_2,\lambda_3)  $ are the mixing parameters and $ \mu=(\mu_1,\mu_2,\mu_3)  $ are the class-specific means. Suppose we have priors

$$\begin{aligned}
\lambda &\sim \text{Dirichlet}(\alpha) \\
\mu_j &\sim \text{Normal}(0,\xi_0) \quad j=1,2,3,
\end{aligned} $$

with $ \alpha=(\alpha_1,\alpha_2,\alpha_3)  $ and where
$ \text{Dirichlet}(\lambda\vert\alpha) = \frac{\Gamma(\sum_j\alpha_j)}{\prod_j\Gamma(\alpha_j)} \prod_j \lambda_j^{\alpha_{j-1}}  $ is the Dirichlet distribution (you can think of it as a generalization of the Beta distribution). The likelihood of the data is

$$  p(y\vert\theta) = \prod_{i=1}^n\Big[\sum_j\lambda_j \text{Normal}(y_i\vert\mu_j,\sigma)\Big]  $$

Assuming that the priors are independent <em>a priori</em>, the posterior is thus

$$\begin{aligned}
p(\theta\vert y) &\propto \prod_{i=1}^n\Big[\sum_j\lambda_j \text{Normal}(y_i\vert\mu_j,\sigma)\Big]\text{Dir}(\lambda\vert\alpha)\text{Normal}(\mu_1\vert0,\xi_0)\text{Normal}(\mu_2\vert0,\xi_0) \\
&\propto \prod_{i=1}^n\left[\sum_j \lambda_j \exp\left(-\frac{1}{2\sigma^2}(y_i-\mu_j)^2\right)\right]\left(\prod_j \lambda_j^{\alpha_j-1}\right)\left[\prod_j \exp\left(-\frac{\mu_j^2}{2\xi_0^2}\right)\right]
\end{aligned} $$

Now, this looks a little bit ugly, especially the sum within the product. But we can introduce a new latent variable $ z_i  $, which indicates class membership, into the model, where $ z_i \sim \text{Categorical}(\lambda)  $, so that $ \Pr[z_i=j]=\lambda_j  $ and $ y_i\vert z_i=j \sim \text{Normal}(\mu_j,\sigma)  $. Including $ z  $ into the set of parameters, we have $ \theta = (\lambda, \mu_1,\mu_2,z)  $, and the posterior will look something like

$$
p(\theta\vert y) \propto \prod_{i=1}^n\left[\lambda_{z_i} \exp\left(-\frac{(y_i-\mu_{z_i})^2}{2\sigma^2}\right)\right]\left(\prod_j \lambda_j^{\alpha_j-1}\right)\left[\prod_j \exp\left(-\frac{\mu_j^2}{2\xi_0^2}\right)\right]
$$

This looks much better. Now, let us consider the conditional distributions. Let $ \bar y_j = \sum_{i:z_i=j}y_i/n_j  $, $ n_j  $ with being the number observations in group $ j=1,2.  $ Dropping unnecessary terms, the conditional distribution of $ \mu_j  $ is simply

$$\begin{aligned}
\mu_j \vert\lambda, \mu_{-j}, z,y &\propto \exp\left(-\frac{1}{2\sigma^2}\sum_{i: z_i=j}\mathbb (y_i-\mu_j)^2\right)\exp\left(-\frac{\mu_j^2}{2\xi_0^2}\right)\\
&= \exp\left(-\frac{1}{2}\left[\sum_{i:z_i=j}\frac{(y_i-\mu_j)^2}{\sigma^2}+\frac{\mu_j^2}{\xi_0^2}\right]\right) \\
&= \exp\left(-\frac{1}{2}\left[\frac{\xi_0^2(\sum_{i:z_i=j}y_i^2-2\mu_jn_j\bar{y}_j+n_j\mu_j^2)+\sigma^2\mu_j^2}{\sigma^2\xi_0^2}\right]\right) \\
&\propto \exp\left(-\frac{1}{2}\left[\frac{\mu_j^2(\sigma^2 + n_j\xi_0^2)-2\mu_jn_j\bar{y}_j\xi_0^2}{\sigma^2\xi_0^2}\right]\right) \\
&\propto \exp\left[-\frac{1}{2}\left(\frac{\sigma^2+ n_j\xi_0^2}{\sigma^2\xi_0^2}\right)\left(\mu_j^2-2\frac{\mu_jn_j\bar{y}_j\xi_0^2}{\sigma^2+n_j\xi_0^2}+ \left(\frac{n_j\bar{y}_j\xi_0^2}{\sigma^2 + n_j\xi_0^2}\right)^2\right)\right] \\
&=\exp\left[-\frac{1}{2}\left(\frac{1}{\xi_0^2}+\frac{n_j}{\sigma^2}\right)\left(\mu_j-\frac{\frac{n_j\bar{y}_j}{\sigma^2}}{\frac{1}{\xi_0^2}+\frac{n_j}{\sigma^2}}\right)^2\right]
\end{aligned}
$$

where $ \mu_{-j}  $ is the vector $ \mu  $ with $ \mu_j  $ excluded. This is the kernel of a normal distribution. As for $ z_i  $, we notice that

$$\begin{aligned}
\Pr[z_i = j\vert\mu,\lambda,y] &\propto \Pr[\mu,\lambda\vert z_i=j,y]\Pr[z_i=j\vert y]\\
&=\text{Normal}(y_i\vert\mu_j,\sigma) \lambda_j
\end{aligned} $$

which is the kernel of a categorical distribution with normalizing constant $ \sum_j \lambda_j \text{Normal}(y_i\vert\mu_j,\sigma)  $. Lastly,

$$\begin{aligned}
\lambda \vert \mu,z,y &\propto \prod_{i}\lambda_{z_i}\prod_j \lambda_j^{\alpha_j-1}
= \prod_j \lambda_j^{n_j}\prod_{j}\lambda_j^{\alpha_j-1}\\
&= \prod_j \lambda_j^{\alpha_i+n_j-1}
\end{aligned} $$

which is the kernel of a Dirichlet distribution.

Thus, we can construct the Gibbs sampler as follows.  At iteration $ t  $,

1\. Sample $ z_i^{t+1}  $ from the discrete distribution on $ \{1,2,3\}  $, where 
    
$$  z_i^{(t+1)} \sim \text{Cat}\left(\frac{\text{Normal}(y_i\vert\mu_1^{(t)},\sigma) \lambda_1^{(t)}}{\sum_j \text{Normal}(y_i\vert\mu_j^{(t)},\sigma) \lambda_j^{(t)} },\frac{\text{Normal}(y_i\vert\mu_2^{(t)},\sigma) \lambda_2^{(t)}}{\sum_j \text{Normal}(y_i\vert\mu_j^{(t)},\sigma) \lambda_j^{(t)} },\frac{\text{Normal}(y_i\vert\mu_3^{(t)},\sigma) \lambda_3^{(t)}}{\sum_j \text{Normal}(y_i\vert\mu_j^{(t)},\sigma) \lambda_j^{(t)} }\right)  $$

2\. Sample 
    
$$  \mu_j^{(t+1)} \sim \text{Normal}\left(\frac{\frac{n_j^{(t)}\bar{y}_j^{(t)}}{\sigma^2}}{\frac{1}{\xi_0^2}+\frac{n_j^{(t)}}{\sigma^2}},\left(\frac{1}{\xi_0^2}+\frac{n_j^{(t)}}{\sigma^2}\right)^{-1/2}\right)  $$

where we note that $ n_j  $ and $ \bar y_j  $ are now random as they depend on $ z  $.

3\. Lastly, update $ \lambda  $ by sampling from 
    
$$  \lambda \sim \text{Dir}(\alpha^{(t)})  $$ 

where $ \alpha_j^{(t)} = \alpha_{0j}+n_j^{(t)}-1  $.

Okay. Let's try it out! First, we load some packages, set the random seed, and generate a dataset.

{% highlight r %}
library('dplyr')
library('data.table')
library('ggplot2')
library('grid')
library('gridExtra')
library('viridis')

# set seed
set.seed(123)

# number of observations
n <- 1000

# parameters
mu <- c(-.3, 0, .3)
sigma <- .15
lambda <- c(.4, .2, .4)

# latent classe sizes
s.class.1 <- round(n * lambda[1])
s.class.2 <- round(n * lambda[2])
s.class.3 <- n - s.class.1 - s.class.2

# latent classes
z <- c(rep(1, s.class.1), rep(2, s.class.2), rep(3, s.class.3))
# outcome
y <- ifelse(z == 1,
        rnorm(s.class.1, mu[1], sigma),
            ifelse(z == 2,
                rnorm(s.class.2, mu[2], sigma),
                rnorm(s.class.3, mu[3], sigma)
            )
    )

df <- data.table(
        y = y,
        z = factor(z, labels = c('Class 1', 'Class 2', 'Class 3')),
        w = sapply(z, function(a) length(y[z==a]) / length(y))
    )
df[, w := w / sum(w)]
{% endhighlight %}

Let's have a look at the generated data.

{% highlight r %}
ggplot(df, aes(x = y, fill = z, weight = w)) +
    geom_histogram(
        alpha = .5, 
        col = 'white', 
        bins = 50, 
        position = 'identity'
    ) +
    scale_fill_viridis(
        name = 'Latent Classes',
        discrete = T,
        begin = 0,
        end = .5,
        option = 'A'
    ) +
    labs(x = 'Data', y = 'Density') +
    theme_bw()
{% endhighlight %}

<img class="alignnone size-full wp-image-1309" src="https://barumpark.files.wordpress.com/2018/03/unnamed-chunk-2-1.png" alt="unnamed-chunk-2-1" width="672" height="480" />

Next, we set up the Gibbs sampler

{% highlight r %}
# number of classes
k<-3
# initial values
xi0 <- 5
alpha0 <- rep(1,k)
mu.t <- rep(0,k)
lambda.t <- rep(1/k, k)
{% endhighlight %}

The code for updating class membership is

{% highlight r %}
update.z <- function() {
    # get probabilities
    p <- sapply(1:k, function(w) dnorm(y,mu.t[w], sigma) * lambda.t[w])
    # draw z
    apply(p, 1, function(w) which(rmultinom(1, 1, w) == 1))
}
{% endhighlight %}

That for the class-means is

{% highlight r %}
update.mu <- function() {

    n.j <- sapply(1:k, function(w) sum(z.t == w))
    j.means <- sapply(1:k, function(w) {
        if (n.j[w] != 0) mean(y[z.t == w]) else 0})
    sapply(1:k, function(w) {
        rnorm(1, 
              (n.j[w] * j.means[w] / sigma^2) / (1 / xi0^2 + n.j[w] / sigma^2),
              (1 / xi0^2 + n.j[w] / sigma^2)^(-1/2))
    })
}
{% endhighlight %}

Lastly, we code the updating function for the mixing parameters:

{% highlight r %}
update.lambda <- function() {
    n.j <- sapply(1:k, function(w) sum(z.t == w))
    MCMCpack::rdirichlet(1, alpha0 + n.j - 1)
}
{% endhighlight %}

Finally, let us set the number of iterations at 15,000 and run the Gibbs sampler!

{% highlight r %}
n.iter =15000
n.burnin = 0
res <- matrix(
    NA, 
    nrow = n.iter - n.burnin, 
    ncol = length(mu.t) + length(lambda.t)
)

# burn-in iterations
for (iter in 1:n.burnin) {

    if (iter %% 2000 == 0)
    message(paste0('Burn-in, Iteration No: ',
            iter, '/', n.burnin))
            
    z.t <- update.z()
    mu.t <- update.mu()
    lambda.t <- update.lambda()
    
}

# sampling
for (iter in 1:(n.iter - n.burnin)) {

    if (iter %% 2000 == 0)
    message(paste0('Sampling, Iteration No: ',
            iter + n.burnin, '/', n.iter))
            
    z.t <- update.z()
    mu.t <- update.mu()
    lambda.t <- update.lambda()

    res[iter, 1:k] <- mu.t
    res[iter, (k+1):(2 * k)] <- lambda.t
}

message('Sampling Done!')

{% endhighlight %}

    Burn-in, Iteration No: 0/0

    Sampling, Iteration No: 2000/15000
    Sampling, Iteration No: 4000/15000
    Sampling, Iteration No: 6000/15000
    Sampling, Iteration No: 8000/15000
    Sampling, Iteration No: 10000/15000
    Sampling, Iteration No: 12000/15000
    Sampling, Iteration No: 14000/15000
    Sampling Done!

A look at the traceplots reveals that the chain reaches a stationary distribution quite fast for this example.

<img class="alignnone size-full wp-image-1315" src="https://barumpark.files.wordpress.com/2018/03/unnamed-chunk-8-11.png" alt="unnamed-chunk-8-1" width="672" height="480" />

{% highlight r %}
df <- as.data.table(
    cbind(1:(n.iter-n.burnin),res)
) %>%
    setnames(
        c('iter',
          paste0('mu[',1:k,']'),
          paste0('lambda[',1:k,']')
          )
    )
    
l.df <- melt(df, id.vars = 'iter')

ggplot(l.df, aes(x = iter, y = value)) +
    geom_line(col = viridis::viridis(1), alpha = .7) +
    theme_bw() +
    facet_wrap(~ variable, ncol = k, scale = 'free') +
    labs(x = 'Iteration', y = 'Sampled Values')

{% endhighlight %}

Let us discard the first 10,000 iterations as burn-in and retain only every second of the remaining half of the iterations (thinning). The posterior distribution of the parameters look like:

{% highlight r %}
df <- df[iter %in% seq(2501, n.iter, 2)]
melt(df, id.vars='iter') %>%
    ggplot(aes(x=value)) +
    geom_histogram(col = viridis::viridis(1),
                   fill = 'white') +
    theme_bw() +
    labs(x = '', y = '')+
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )+
    facet_wrap(~ variable, ncol = k, scale = 'free')
{% endhighlight %}

<img class="alignnone size-full wp-image-1311" src="https://barumpark.files.wordpress.com/2018/03/unnamed-chunk-9-1.png" alt="unnamed-chunk-9-1" width="672" height="480" />

So, it seems that the code is working! In this example, we have not encountered what is often referred to as <em>label switching</em>. Namely, we could label the classes as class 1, class 2, and class 3, as we did here, but we could also label them in the order 3, 2, 1 and switch the labels of the corresponding means and mixing parameters in the same way. Notice that such relabeling will lead to the exact same likelihood and since we are assigning the same priors for the parameters associated with each label, the posterior probabilities for each of these configurations will be same as well. In fact, if we have $K$ components, we know _a priori_ that the posterior will have at least $K!$ modes that are equally likely, one mode for each permutation of the labels. That our MCMC chains does not show a switching between any of the modes is, therefore, an indication that the sampler is not exploring the entire posterior. On the other hand, if the chains switch between labels, we need additional methods to relabel the classes. Otherwise, we would have to give up any inference regarding label-specific parameters, such as the means in this example. There is in fact a growing literature on label switching in Bayesian mixture models and algorithms have been developed to ensure that the the modes corresponding to all permutations of the labels are sufficiently explored.
