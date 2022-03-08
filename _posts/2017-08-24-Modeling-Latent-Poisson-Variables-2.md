---
layout: post
title: Modeling Latent Poisson Variable (pt. 2)
date: 2017-08-23 21:45
author: baruuum
comments: true
categories: [Bayesian, Quant Stuff, STAN]
---

**DISCLAIMER**: _This is a blog post from my graduate school years. The blog is no longer maintained and the material might include typos._

This is a continuation of the [previous post](/blog/2017/Modeling-Latent-Poisson-Variables) about estimating the rate parameter of a Poiosson random variable $Y^* $, which is only partially observed. This time we'll try to estimate $\lambda$ using Bayesian estimation in `STAN` via the `rstan` package in `R`. As before, let us fist consider a situation in which we can observe the latent variable $Y^* \sim \text{Poisson}(\lambda)$. In this case the estimation is simple.

But first, let us generate exactly the same data as in the previous post.

{% highlight r %}
set.seed(321)
n = 5000
lambda = 2
y_star = rpois(n, lambda)
y = ifelse(y_star == 0, 1,
         ifelse(y_star == 1, 2,
                ifelse(1 < y_star & y_star < 5, 3, 4)))
{% endhighlight %}

Next, we load the `rstan` package

{% highlight r %}
library('rstan')
{% endhighlight %}

and define the model

{% highlight r %}
model.string <- '
data {
   int N;
   int y[N];
}

parameters {
   real omega;
}

transformed parameters {
   real lambda;

   lambda = 1 / omega;

}

model {

   omega ~ uniform(0, 1);  // prior
   y ~ poisson(lambda); // likelihood

}
'
{% endhighlight %}

The model is straightforward. There are `N` observations and the only data imput is a length-`N` vector of responses, `y`. The responses are assumed to be iid samples from a Poisson distribution with parameter $\lambda$. As $\lambda \in (0,\infty)$, we use $1/\omega$ as the prior distribution for $\lambda$, where $\omega\sim \text{Uniform}(0,1)$. (Edit: some reflection makes me think that this is not a good prior distribution as it will assign zero density to the range of values between zero and one.) Next, we fit the model:

{% highlight r %}
# data
stan.dat <- list(N=length(y_star), y = y_star)

# start sampling
fit <- stan(model_code = model.string,
            data = stan.dat,
            iter = 5000,
            chains = 1,
            warmup = 500,
            refresh = 1000,
            thin = 2)
{% endhighlight %}

    ## SAMPLING FOR MODEL 'a7c53550e0ae4b65ea895e70db308792' NOW (CHAIN 1).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 5000 [  0%]  (Warmup)
    ## Iteration:  501 / 5000 [ 10%]  (Sampling)
    ## Iteration: 1500 / 5000 [ 30%]  (Sampling)
    ## Iteration: 2500 / 5000 [ 50%]  (Sampling)
    ## Iteration: 3500 / 5000 [ 70%]  (Sampling)
    ## Iteration: 4500 / 5000 [ 90%]  (Sampling)
    ## Iteration: 5000 / 5000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 0.633 seconds (Warm-up)
    ##                5.154 seconds (Sampling)
    ##                5.787 seconds (Total)


We can look into a summary of the posterior distribution by using the `summary` command as follows:

{% highlight r %}
summary(fit)$summary
{% endhighlight %}

    ##                 mean      se_mean          sd 
    ## omega      0.4976819 0.0001336453 0.004981163 
    ## lambda     2.0095168 0.0005396766 0.020117482 
    ## lp__   -3036.3338851 0.0172031543 0.681428219 
    ##                  50%     n_eff      Rhat
    ## omega      0.4976294  1389.165 1.0008747
    ## lambda     2.0095275  1389.569 1.0008913
    ## lp__   -3036.0819005  1569.004 0.9998323

The posterior mean and median are, again, close to $\lambda = 2$. But notice that here we have assumed that $Y^* $ is observed. If not, how can we estimate $\lambda$ with the categorized data? A little bit of thought reveals that this is not a straightforward problem. The main problem consists of expressing the likelihood

$$ L(\lambda ;  \mathbf{y}) = p(\mathbf y\vert \lambda) = \prod_{i=1}^n \left[ \sum_{l=1}^K \sum_{y_i^{*} \in I_l} \frac{\lambda^{y_i^{*}}}{y_i^{*}!}\exp(-\lambda)\times \mathbb I_{I_l}(y_i^{*})\right]. $$

As this is no standard PMF, there exists no standard command such as `y ~ poisson(lambda)`, which we have used above. Rather we have to calculate the target log-density directly. By the target density, we mean the joint density of the parameters from which the MCMC sampler samples. `STAN` has an built-in object called `target` which stands for the "log" target density, $\ln p(\theta\vert  \mathbf y)$, where $\theta$ is the parameter vector of interest and $\mathbf y$ is the data vector. To fix ideas, recall that the posterior is proportional to the likelihood times the prior. Thus

$ \ln p(\theta\vert  \mathbf y) = \ln K + \ln p(\theta) + \ln p(\mathbf y \vert \theta), $

where $K$ is a constant that does not depend on the parameters. As the posterior is a linear function of the log-likelihood and the log-prior density, we can add the contribution of each observation to the log-likelihood one-by-one in order to obtain the posterior density we desire. For example,

    y ~ normal(mu, sigma)

is equivalent to

    for (i in 1:N) y[i] ~ normal(mu, sigma)

and does the same thing in `STAN` as

    target += normal_lpdf(y | mu, sigma)

which is, again, equivalent to

    for (i in 1:N) target += normal_lpdf(y[i] | mu, sigma).

The meaning of the first statement is quite clear, we are asserting that the likelihood is the pdf of a $N$ iid $\text{Normal}(\mu, \sigma)$ variables. With this statement `STAN` will calculate the log-likelihood for us, which is then used to sample from the posterior. But, rather than writing it in a "sampling statement" as in the first line, we can simply take the log of the density function of a $\text{Normal}(\mu, \sigma)$ variable, evaluate it at each data point, and sum these values up. This is what we are doing in the second and the third line. `normal_lpdf(y | mu,sigma)` is the log of the normal density with parameters `mu` and `sigma` evaluated at `y`. By putting this expression on the RHS of the syntax, we require `STAN` to increment the log target density by this amount, which is equivalent to the adding each observation's contribution to the log-likelihood.

Thus, the model can be written as

{% highlight r %}
model.string.cat <- '
data {
   int N;
   int y[N];
}

parameters {
   real omega;
}

transformed parameters {
   real lambda;

   lambda = 1 / omega;

}
model {

   omega ~ uniform(0, 1);  // prior
   for (nn in 1:N) { // likelihood
      if (y[nn] == 1) {
         target += poisson_lpmf(0 | lambda);
      } else if (y[nn] == 2) {
         target += poisson_lpmf(1 | lambda);
      } else if (y[nn] == 3) {
         target += log(exp(poisson_lcdf(4 |lambda))
                       - exp(poisson_lcdf(1| lambda)));
      } else {
         target += poisson_lccdf(4 | lambda);
      }
   }

}
'
{% endhighlight %}

Here `poisson_lpmf`, `poisson_lcdf`, and `poisson_lccdf` are respectively the log Poisson PMF, the log CDF, and the log complementary CDF. That is if $f(x)$ is the PMF and $F(x)$ the CDF of a Poisson distribution evaluated at $x$, then `poisson_lpmf` is $\ln f(x)$, `poisson_lcdf` is $\ln F(x)$, and `poisson_lccdf` is $\ln[1 − F(x)]$.

Fitting the model to the data, we get the following results:

{% highlight r %}
#data
stan.dat.2 <- list(N = length(y), y = y)

# start sampling
fit <- stan(model_code = model.string.cat,
            data = stan.dat.2,
            iter = 5000,
            chains = 1,
            warmup = 500,
            refresh = 1000,
            thin = 2)
{% endhighlight %}

    ## SAMPLING FOR MODEL 'f65ed971427f4b91bba4006bf17bc125' NOW (CHAIN 1).
    ## 
    ## Gradient evaluation took 0.003 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 30 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 5000 [  0%]  (Warmup)
    ## Iteration:  501 / 5000 [ 10%]  (Sampling)
    ## Iteration: 1500 / 5000 [ 30%]  (Sampling)
    ## Iteration: 2500 / 5000 [ 50%]  (Sampling)
    ## Iteration: 3500 / 5000 [ 70%]  (Sampling)
    ## Iteration: 4500 / 5000 [ 90%]  (Sampling)
    ## Iteration: 5000 / 5000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 9.301 seconds (Warm-up)
    ##                91.2 seconds (Sampling)
    ##                100.501 seconds (Total)

{% highlight r %}
summary(fit)$summary
{% endhighlight %}

    ##                 mean      se_mean         sd 
    ## omega      0.4995896 0.0001350554 0.00538462 
    ## lambda     2.0018757 0.0005415327 0.02159037 
    ## lp__   -5601.2570369 0.0181856058 0.67236024 
    ##                  50%    n_eff      Rhat
    ## omega      0.4997208 1589.594 0.9996615
    ## lambda     2.0011176 1589.537 0.9996567
    ## lp__   -5600.9989037 1366.937 0.9997972


Again, extremely close!
