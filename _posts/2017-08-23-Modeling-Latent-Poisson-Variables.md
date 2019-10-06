---
layout: post
title: MLE of Latent Poisson Variable
date: 2017-08-23 05:48
author: baruuum
comments: false
categories: [MLE, Quant Stuff]
---

Suppose we are interested in a random variable $$Y^*$$ which follows a Poiosson distribution with parameter $$\lambda$$, i.e.,

$$ Y^* \sim \text{Poisson}({\lambda}).$$

Obviously, we are not able to observed $$\lambda$$. To make things even worse, suppose that can't observe $$Y^*$$ either, but only a categorized form of $$ Y^* $$, which we denote by $$Y$$, defined as $$ Y = 1 $$ if $$ Y^* = 0$$, $$Y = 2$$ if $$Y^* = 1$$, $$ Y = 3 $$ if $$ Y^* \in \{2,3,4\} $$, and $$ Y = 4 $$ if $$ Y^* > 5 $$.

We have a data vector $$ \mathbf{y} = [y_1,y_2,...,y_n], y_i \in \{1,2,3,4\}, $$ of length $$n$$. Now, how to estimate $$\lambda$$?

If we were able to observe $$Y^*$$, then the likelihood can be easily formed. Given that the observations are independently drawn and are coming from the same distribution, we have

$$ L(\lambda; \mathbf{y}^*) = \prod_{i=1}^n \frac{\lambda^{y_i^*}}{y_i^*!}\exp(-\lambda). $$

Let us translate this into the likelihood of $$\lambda$$ given $$ \mathbf y $$. $$Y^* $$ can fall into any of the four intervals defined above. The probability that $$Y^* $$ falls in to the first category ($$Y^* = 0$$) is $$\exp(−\lambda)$$ (simply plug in 0 into the Poisson PMF). Similarly, the probability of $$Y^*$$ falling into the third category is $$\sum_{k=2}^4 (\lambda^k/k!)\exp(-\lambda)$$. In general, therefore, we can define the likelihood of $$\lambda$$ as

$$
L(\lambda ; \mathbf{y}) = \prod_{i=1}^n \left[ \sum_{l=1}^K \sum_{y_i^*\in I_l} \frac{\lambda^{y_i^*}}{y_i^*!}\exp(-\lambda)\times \mathbb I_{I_l}(y_i^*)\right]
$$

where $$I_l$$ is the $$l$$th interval into which $$Y^* $$ can fall, $$K$$ is the total number of intervals, and $$ \mathbb I_A(x) $$ is an indicator function which is 1 if $$ x\in A $$ and zero otherwise. The likelihood looks complicated, but is actually straightforward to interpret. We work from the inside outwards. First, in the middle of the equation is the Poisson PMF. We sum the PMF (which gives the probability of observing $$y_i^* $$) for all $$y_i^* $$ values which are in the $$l$$th interval. The outer sum is over all such intervals; in this example we have 4 of them. Thus, without the indicator function, the term in the square brackets will be always one. What the indicator function does is singling out those values that are in the interval $$I_l$$.

But did we do everything right...? If the usual regularity assumptions are satisfied, which is the case in this example, the MLE of $$ \lambda $$ should be quite close to the true value $$ \lambda $$ if the sample size is reasonably large. So, let us simulate a data for which we know the exact data-generating process and, thus, the "true" values of the parameters.

{% highlight r %}
# set seed
set.seed(321) 
# no. of obs.
n = 5000
# par. value
lambda = 2
# generate latent var.
y_star = rpois(n, lambda)
# generate observed var.
y=ifelse(y_star == 0, 1,
         ifelse(y_star == 1, 2,
                ifelse(1 < y_star & y_star < 5, 3, 4)))
{% endhighlight %}

As the log function is strictly increasing on $$ \mathbb R_+ $$, maximizing the likelihood is equivalent to maximizing the log-likelihood:

$$
\ln L(\lambda ; \mathbf y) = \sum_{i=1}^n \ln \left[\sum_{l=1}^K\sum_{y_i^* \in I_l} \frac{\lambda^{y_i^*}}{y_i^*!}\exp(-\lambda)\times \mathbb I_{I_l}(y_i^*\in I_l)\right].
$$

So, let's write a function

{% highlight r %}
log.lik = function(lambda,dat=y) {
   # empty vector to store results
   ll.vec = numeric(n)
   for (ii in 1:n) { # generate likelihood
      if (y[ii] == 1) {
         ll.vec[ii] = dpois(0, lambda)
      } else if (y[ii] == 2) {
         ll.vec[ii] = dpois(1, lambda)
      } else if (y[ii] == 3) {
         ll.vec[ii] = sum(dpois(2:4, lambda))
      } else {
         ll.vec[ii] = 1 - ppois(4, lambda)
      }
   }
   # log the results and then sum
   return(sum(log(ll.vec)))
}
{% endhighlight %}

which log-likelihood function. To maximize this function numerically we use the 'BFGS' algorithm offered in the `optim` function:

{% highlight r %}
init = runif(1, 0, 5) # reasonable initial value
res <- optim(par = init,
             fn = log.lik,
             dat = y, 
             method = 'L-BFGS-B',
             lower = 0,
             control = list(fnscale = -1))

cat(paste0('Convergece Code: ', res$convergence))
{% endhighlight %}

     # Convergece Code: 0


A convergence code of 0 implies that the algorithm has successfully converged to a local maximum. Our Maximum Likelihood estimate is:

{% highlight r %}
cat(paste0('Parameter Estimate: ', round(res$par,3)))
cat(paste0('Log-Likelihood: ', round(res$value,3)))
{% endhighlight %}


     # Parameter Estimate: 2.002
     # Log-Likelihood: -5599.384


That is indeed close to the true parameter value of $$ \lambda=2 $$. Thus, it seems that it works!
