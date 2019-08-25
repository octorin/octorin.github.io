---
layout: post
title: Misleading Textbook Example
date: 2018-08-16 03:20
author: baruuum
comments: true
categories: [Bayesian, Quant Stuff]
---
One of the favorite examples with which many Bayesian introductory textbooks start is updating the belief of the probability a coin will land up heads. The example goes something like this. Let $ \theta  $ be the probability that a coin will land up heads. We have some prior beliefs regarding $ \theta  $. As a probability has to lie between zero and one, a natural choice is to express our beliefs with a Beta distribution, which has probability density function

$$  f(x\vert\alpha,\beta) = \frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\Gamma(\beta)}x^{\alpha - 1}x^{\beta-1}, \qquad x\in [0,1], \alpha>0,\beta>0,   $$

where $ \Gamma(\cdot)  $ is the Gamma function, defined by the improper integral $ \Gamma(z) = \int_0^\infty u^{z-1}e^{-u} du  $ which converges absolutely for real $ z  $. The mean of the beta distribution is

$$ \begin{aligned}
E[X] &= \int_0^1 x f(x\vert\alpha,\beta) dx  = \int_0^1 \frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\Gamma(\beta)}x^{\alpha}x^{\beta-1}dx \\
&=  \frac{\Gamma(\alpha+1)\Gamma(\alpha + \beta)}{\Gamma(\alpha+\beta+1)\Gamma(\alpha)}\int_0^1 \frac{\Gamma(\alpha+\beta+1)}{\Gamma(\alpha+1)\Gamma(\beta)}x^{\alpha}x^{\beta-1}dx \\
&= \frac{\alpha}{\alpha + \beta},
\end{aligned}  $$

using the identity $ \Gamma(u+1) = u\Gamma(u)  $ and noting that the integrand in the second line is just the $ \text{Beta}(\alpha+1,\beta)  $ pdf which has to integrate to one. So if we choose $ \theta \sim \text{Beta}(\alpha,\beta)  $ with $ \alpha=\beta  $, the prior mean will be $ 1/2  $.

Now, suppose our prior beliefs can be expressed by the $ \text{Beta}(4,4)  $ distribution. The $ \text{Beta}(4,4)  $ pdf looks like this:

{% highlight r %}
# grid
x <- seq(0, 1, .0005)

#plot
plot(x, dbeta(x, 4, 4), 
     ylab = '', 
     xlab = expression(theta), 
     col = viridis::viridis(1, begin = .5, option = 'A'), 
     type = 'l',
     axes = F,
     bty = 'n')
axis(1, at = seq(0, 1, .25))
title('Beta(4,4) Prior Distribution')
{% endhighlight %}

<center>
<img src="/assets/img/misleading1.png" width="500" height="400" />
</center>

Next, supose we observe $ n=10  $ independent coin tosses, out of which 2 land heads. Assuming that the probability of heads is constant across the trials, this implies that $ y \vert \theta \sim \text{Binomial}(n, \theta)  $. Thus, the likelihood of the data is

$$   p(y\vert\theta) = \binom{10}{2} \theta^{\sum_{i=1}^{10} y_i}(1-\theta)^{10-\sum_{i=1}^{10} y_i},   $$

with $ \sum_{i=1}^{10} y_i = 2.  $ By Bayes' Rule, the posterior distribution is proportional to the likelihood times the prior:

$$  \begin{aligned}
p(\theta\vert y) &\propto  \underbrace{\binom{10}{2}\theta^{\sum_{i=1}^{10} y_i}(1-\theta)^{10-\sum_{i=1}^{10} y_i}}_{\text{Likelihood}} \times \underbrace{\frac{\Gamma(8)}{\Gamma(4)\Gamma(4)}\theta^3 (1-\theta)^3}_{\text{Prior}} \\
&\propto \theta^{2}(1-\theta)^{10-2}\theta^3 (1-\theta)^3\\
&=\theta^{5}(1-\theta)^{11}.
\end{aligned}  $$

As the posterior distribution has to be a <em>distribution</em>, which integrates to one, all that is left is figuring out the proportionality constant $ K  $ that satisfies

$$  K\int_{[0,1]}\theta^{5}(1-\theta)^{11}d\theta = 1.   $$

But a close look into the integrand reveals that it looks much like a Beta pdf, except that the constant with the Gamma functions is missing. Thus, $ K = \Gamma(6+12)/[\Gamma(6)\Gamma(12)]  $ will satisfy the equation, showing that the posterior distribution is again a Beta distribution with shape parametres $ \alpha_{\text{post}} = 6  $ and $ \beta_{\text{post}} = 12  $.

In general, when we have $ n  $ independent Bernoulli trials out of which $ k  $ are a "success" and a $ \text{Beta}(\alpha,\beta)  $ prior, the posterior distribution will be a Beta distribution with shape parameters $ \alpha_{\text{post}} = \alpha + k  $ and $ \beta_{\text{post}} = \beta + (n - k)  $, as can be read out from the derivation above. Thus, for fixed $ n  $, the more successes we have the higher will be the posterior mean, which makes sense.

Notice that both the prior and the posterior distribution belong to the same family of distributions. The prior that satisfied this condition is called the <em>conjugate prior</em> for the likelihood. Thus, the Beta distribution is the conjugate prior for the Binomial likelihood.

The prior, likelihood, and posterior from the example above look something like this:

{% highlight r %}
# colors
col.3 <- viridis::viridis(3, end = .8, option = 'A')

# prior
plot(x, dbeta(x, 4, 4), 
     ylab = '', 
     xlab = expression(theta), 
     col = col.3[1], 
     type = 'l',
     axes = F,
     ylim = c(0, 4),
     bty = 'n')
axis(1, at = seq(0,1,.25))

# likelihood
curve(dbinom(2, 10, x), add = T, col = col.3[2])

# posterior
curve(dbeta(x, 6, 12), add = T, col = col.3[3])
legend('topright', col=col.3, lty = rep(1,3),
       legend = c('Prior','Likelihood','Posterior'),
       bty='n')
title('Prior, Likelihood, and Posterior')
{% endhighlight %}

<center>
<img src="/assets/img/misleading2.png" width="500" height="400" />
</center>

Recall that <strong>the likelihood is not a probability distribution of $ \theta  $: it is the Binomial pmf treated as a function of $ \theta  $.</strong> It will sum to one if we sum of the number of successes from $ 0  $ to $ n=10  $ with fixed $ \theta  $; but it will **not** integrate to one if we integrate the likelihood function with respect to $ \theta  $ over $ (0,1)  $ for fixed number of successes and trials. Thus, the area under the purple curve in the figure is not equal to one. For the example case we can calculate the area under the curve explicitly as 

$$ \binom{10}{2} \int_{(0,1)} \theta^2(1-\theta)^8 d\theta = \frac{1}{11}.$$

So, what people often do is to scale the likelihood such that the area under the curve becomes equal to one. Using this _trick_, we get the well-known plot

{% highlight r %}
# prior
plot(x, dbeta(x, 4,4), 
     ylab = '', 
     xlab = expression(theta), 
     col = col.3[1], 
     type = 'l',
     axes = F,
     ylim = c(0, 4),
     bty = 'n')
axis(1, at=seq(0, 1, .25))

# scaled likelihood
curve(dbeta(x, 3, 9), add = T, col = col.3[2])

# posterior
curve(dbeta(x, 6, 12), add = T, col = col.3[3])
legend('topright', col = col.3, lty = rep(1,3),
       legend = c('Prior', 'Likelihood', 'Posterior'),
       bty = 'n')
title('Prior, (Scaled) Likelihood, and Posterior', cex.main = .9)
{% endhighlight %}

<center>
<img src="/assets/img/misleading3.png" width="500" height="400" />
</center>

I think that this plot is quite misleading, since the curve representing the likelihood function is in fact representing something else. So, it creates the false impression that $ p(y\vert\theta)  $ is a density function of $ \theta  $ and not $ y$.

A second note, although not a "misleading" one, is that most of the time in sociological research, you'll never observe such a plot. This is because we have often $ n \gg 500  $. Even with $ n= 500  $, the likelihood will dominate the posterior, so that the influence of the prior on the posterior will be minimal: the prior will provide only 1/500 of the information that the likelihood provides.

To demonstrate this, consider a set of prior densities with different means. We might generate such a set of priors by reparameterizing the Beta distribution with a mean parameter $ \mu = \alpha/(\alpha + \beta)  $ and a "precision" parameter $ \nu = \alpha + \beta  $. Then, $ \alpha = \mu\nu  $ and $ \beta = (1-\mu)\nu  $.

{% highlight r %}
# function to transform mean and precision to shape params
mean.to.shape <- function(mu, nu) {
  c(alpha = mu * nu, beta = (1 - mu) * nu)
}

# set plotting parameters
par(mfrow = c(1, 2), mar = c(4, .1, 4, .5))

# sample size 
n <- 500

# proportion of "successes"
p <- .2

# plot of priors
plot(x, dbeta(x, 4, 4), 
     ylab = '', 
     xlab = expression(theta[prior]), 
     col = 'blue', 
     type = 'n',
     axes = F,
     ylim = c(0, 10),
     bty = 'n')
axis(1, at = seq(0, 1, .2))

# precision
nu <-  15
# means
means <- seq(.1, .9, .1)
# colors
cols <- viridis::viridis(
  length(means), begin = 0, end = .8, option = 'A'
)
# counter
ii <-  0
for (mu in means) {
  ii <- ii + 1
  # get shape parameters
  shape.pars <- mean.to.shape(mu, nu)
  alpha <- shape.pars[1]
  beta <- shape.pars[2]
  # plot prior
  curve(dbeta(x, alpha, beta), 
        add = T, 
        col = cols[ii])
  # add mode of prior distribution
  if (alpha > 1 & beta <1 ) {
    abline(v = 1, col = cols[ii], lty = 2)
  } else {
    abline(v = (alpha - 1) / (alpha + beta - 2),
           col = cols[ii],
           lty = 2)
  }
}
title('Set of (Strong) Priors')

# plot of priors
plot(x, dbeta(x, 4, 4), 
     ylab = '', 
     xlab = expression(theta[posterior]), 
     col = 'blue', 
     type = 'n',
     axes = F,
     ylim = c(0, 25),
     bty = 'n')
axis(1, at = seq(0, 1, .2))

# counter
ii <-  0
for (mu in means) {
  ii <- ii + 1
  shape.pars <- mean.to.shape(mu, nu)
  alpha <- n * p + shape.pars[1]
  beta <- n * (1 - p) + shape.pars[2]
  curve(dbeta(x, alpha, beta),
        add = T, col = cols[ii])
  abline(v = (alpha - 1) / (alpha + beta - 2),
         col = cols[ii],
         lty = 2)
}
title(main = 'Corresponding Posteriors')
mtext(paste0('(n = ', n, ')'))
{% endhighlight %}

<center>
<img src="/assets/img/misleading4.png" width="800" height="300"/>
</center>

Notice that most of these prior distributions would be considered very <em>informative</em>. Indeed, starting from the left, the prior probability that $ \theta  $ lies at $ .2  $ (the MLE) or below is, respectively,

{% highlight r %}
ii <-  0
probs <- numeric(length(means))
for (mu in means) {
  ii <- ii + 1
  shape.pars <- mean.to.shape(mu, nu)
  alpha <- shape.pars[1]
  beta <- shape.pars[2]
  probs[ii] <- pnorm(.2, alpha, beta)
}
print(probs, digits = 4)
{% endhighlight %}

    ## [1] 4.616e-01 4.078e-01 3.411e-01 2.596e-01 1.652e-01 7.123e-02 1.104e-02
    ## [8] 4.189e-05 3.768e-19

The posterior distributions, on the other hand, look quite similar. We can calculate also some posterior quantiles:

{% highlight r %}
# counter
ii <-  0
post.intervals <- matrix(NA, nrow = length(means), ncol = 3)
for (mu in means) {
  ii <- ii + 1
  shape.pars <- mean.to.shape(mu, nu)
  alpha <- n * p + shape.pars[1]
  beta <- n * (1 - p) + shape.pars[2]
  post.intervals[ii, ] <- qbeta(c(.025, .5, .975), alpha, beta)
}
colnames(post.intervals) <- c('q025', 'q50', 'q975')
print(post.intervals, digits = 4)
{% endhighlight %}

    ##         q025    q50   q975
    ##  [1,] 0.1639 0.1967 0.2325
    ##  [2,] 0.1666 0.1996 0.2356
    ##  [3,] 0.1693 0.2025 0.2387
    ##  [4,] 0.1720 0.2054 0.2418
    ##  [5,] 0.1748 0.2084 0.2448
    ##  [6,] 0.1775 0.2113 0.2479
    ##  [7,] 0.1802 0.2142 0.2510
    ##  [8,] 0.1829 0.2171 0.2541
    ##  [9,] 0.1857 0.2200 0.2572

As we can see, the posterior median <em>is</em> moved towards the prior median (or mode). However, even with this set of strong priors, the smallest and largest posterior median differ only by about $ .02  $ and the 95% posterior intervals are all overlapping with one another. Of course, textbooks also explain that the likelihood will dominate when $ n  $ is large. But, given the small differences in the posterior, it is hardly understandable why "applied" sociologists would be skeptical of weakly informative priors, as even strong priors lead to similar conclusions with moderately large sample sizes.

Thus, it is true that Bayesians with different priors will disagree (in a quite technical sense) with one another, even if they have nearly infinite data. Yet, it is probably true as well that the difference will not be very meaningful in applied settings if the sample is large and weakly informative priors are used. An exception would be the case in which the priors assign zero probability to a region within the support of the parameter. But if so, these priors could be hardly called "weakly" informative.

Of course, given these results, the question arises why using Bayesian methods at all or whether Bayesian approaches have no advantages over Frequentists approaches. There is no straightforward answer, but I think that it is precisely where the prior is informative that Bayesian methods are particularly useful. We can think of the prior as providing us with an additional piece of information that can regularize our inference. In the case of parameters that are estimated with abundant data, the influence of this additional piece of information will be negligible, as shown in the example above. But for parameters for which we have only little information, the prior will shrink the estimates towards more reasonable values, assuming that our prior is centered at a reasonable value. 

Think of a multilevel modeling context with individuals nested within groups. For groups with only one or two individuals, we would have very little information to estimate the group-level mean. If the group is heterogeneous, our estimate can land quite everywhere simply due to sampling variability. Yet, by assigning a distribution over the group-level means (i.e., a prior), we can shrink the estimate towards the grand mean, which will prevent out-of-the-world estimates. This is the mechanism through which multilevel models, in general, out-perform fixed-effects models in prediction tasks. It comes from the shrinkage induced by our prior. And we'll have more shrinkage in those areas of the parameter space where information is lacking. 