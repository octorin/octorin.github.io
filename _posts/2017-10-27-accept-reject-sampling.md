---
layout: post
title: Accept-Reject Sampling
date: 2017-10-27 02:32
author: baruuum
comments: false
categories: [Accept-Reject Sampling, Monte Carlo Methods, Quant Stuff, Sampling, von Neumann]
---
Suppose we want to sample from a density {% katex %} p(\theta)  {% endkatex %}, where {% katex %} \theta \in \Theta {% endkatex %} and {% katex %} \Theta {% endkatex %} is the parameter space (here assumed to be unidimensional). Suppose further that {% katex %} p  {% endkatex %} has two properties

<ul>
<li>We can evalulate {% katex %} p  {% endkatex %} for every {% katex %} \theta \in \Theta  {% endkatex %}.</li>
<li>But we are not able to sample directly from {% katex %} p  {% endkatex %}. </li>
</ul>

If we are able to invert the the cdf of {% katex %} \theta  {% endkatex %}, we might use the inverse-probability transform to generate random samples from {% katex %} p  {% endkatex %}. But, suppose that this is not possible. What should we do?

The accept-reject sampling method can be used in such situations if we can find a function {% katex %} g  {% endkatex %} such that 

1. {% katex %} g \ge p {% endkatex %} on {% katex %} \Theta  {% endkatex %} and 
2. {% katex %} g(\theta) = k \cdot h(\theta)  {% endkatex %}, where {% katex %} k  {% endkatex %} is a constant and {% katex %} h  {% endkatex %} is a probability distribution from which we can sample directly. 

Note that {% katex %} g  {% endkatex %} will be in general not a probability distribution as it will not integrate to 1 but instead to {% katex %} k  {% endkatex %}. Given that we have found such a {% katex %} g  {% endkatex %}, the algorithm works as follows:

For a large number {% katex %} N  {% endkatex %},

<ol>
<li>sample {% katex %} x \sim h(\theta)  {% endkatex %}</li>
<li>sample {% katex %} u \sim \text{Unif}(0,1)  {% endkatex %},</li>
<li>calculate {% katex %} r = p(x)/g(x) = p(x)/(k\cdot h(x))  {% endkatex %}</li>
<li>if {% katex %} u \le r  {% endkatex %}, accept the draw and set {% katex %} \tilde{\theta}^{(n)}= x  {% endkatex %} </li>
<li>else, return to step 1. </li>
</ol>

Why does the algorithm work? We will first lay out the intuition and thereafter a more rigorous statement under the assumption that the pdf or pmf of the target distribution exists. Actually, the intuition behind is quite simple: the probability that a proposal {% katex %} x  {% endkatex %} is chosen is {% katex %} h(x)  {% endkatex %}. Given that we have {% katex %} x  {% endkatex %}, the probability that it will be accepted is {% katex %} r=p(x)/(kh(x))  {% endkatex %}. Multiplying both gives us the probability that {% katex %} x  {% endkatex %} is chosen <em>and</em> accepted, which is {% katex %} p(x)/k\propto p(x)  {% endkatex %}. Thus, by repeating the steps in the algorithm, we are sampling proportional to the target density.

Let us express this intuition more formally. Consider the probability {% katex %} \Pr[\tilde{\theta}\le q]  {% endkatex %} for some constant {% katex %} q  {% endkatex %} and note that {% katex %} X  {% endkatex %} and {% katex %} U  {% endkatex %} are independent, so that their joint density is simply {% katex %} f(x,u)=h(x)  {% endkatex %}. Further, this probability can be written as

{% katex display %} \begin{aligned} 
\Pr[\tilde{\theta} \le q ] &= \Pr[x \le q | u \le r] = \frac{\Pr[x \le q \text{ and } u \le r]}{\Pr[u \le r]} \\
&= \frac{\Pr[x \le q \text{ and } u \le p(x)/(kh(x))]}{\Pr[u \le p(x)/(kh(x))]} \\
&= \frac{\int_{-\infty}^q \int_{0}^{p(x)/kh(x)} h(x) dudx}{\int_{-\infty}^{\infty}\int_{0}^{p(x)/kh(x)}h(x)dudx }
\end{aligned}{% endkatex %}

But as {% katex %} u \sim \text{Unif}(0,1)  {% endkatex %}, we have that {% katex %} \int_0^{p(x)/kh(x)} du = \frac{p(x)}{k h(x)}  {% endkatex %}. Thus,

{% katex display %}   \frac{\int_{-\infty}^q \int_{0}^{p(x)/kh(x)}du h(x) dx}{\int_{-\infty}^{\infty}\int_{0}^{p(x)/kh(x)}duh(x)dx} = \frac{\frac{1}{k}\int_{-\infty}^q p(x) dx}{\frac{1}{k}\int_{-\infty}^\infty p(x)dx }=\int_{-\infty}^q p(x)dx.    {% endkatex %}

Therefore, we have {% katex %} \Pr[\tilde \theta \le q ] = \int_{-\infty}^q p(\theta)d\theta  {% endkatex %} as desired.

Note that {% katex %} p(\cdot)  {% endkatex %} appears both in the numerator and the denominator in the derivation. This implies that the density {% katex %} p(\cdot)  {% endkatex %} needs to be known only up to a proportionality constant. The important requirement for the algorithm is that we find a density {% katex %} g = k\cdot h  {% endkatex %} that dominates {% katex %} p  {% endkatex %} on {% katex %} \Theta  {% endkatex %}.

Let us try the algorithm out. Suppose we want to sample from a chi-squared distribution with 2 degrees of freedom. Of course, we know that the distribution is nothing but the distribution of the sum of two squared standard Normal variables, but let us assume we do not know how to sample from it directly. We might use an exponential density with rate parameter {% katex %} \lambda=1/2  {% endkatex %}, {% katex %} h(\theta; \lambda) = \lambda \exp(-\lambda \theta)  {% endkatex %}, and choose the scale constant {% katex %} k=2  {% endkatex %}, which will dominate a chi-squared density with 2 degrees of freedom, as shown in the following graph:

{% highlight r linenos%}
library(data.table)
library(ggplot2)
library(gridExtra)

df <- data.table(theta=seq(0,15,.1))
df[,chi2:=dchisq(theta,2)]
df[,g.x:=exp(-.5*theta)]
df[,l.chi2:=log(chi2)]
df[,l.g.x:=log(g.x)]

cols <- viridis::viridis(2, end=.5)

g1 <- ggplot(df, aes(x = theta)) +
    geom_line(aes(y = chi2), col = cols[1]) +
    geom_line(aes(y = g.x), col = cols[2]) +
    theme_bw() +
    labs(x = expression(theta),
         y = expression(f(theta)))

g2 <- ggplot(df, aes(x = theta)) +
    geom_line(aes(y = l.chi2), col = cols[1]) +
    geom_line(aes(y = l.g.x), col = cols[2]) +
    theme_bw() +
    labs(x = expression(theta),
         y = expression(log(f(theta))))

grid.arrange(g1, g2, nrow = 1)
{% endhighlight %}


<img src="{{ site.baseurl }}/assets/img/arsampling1.jpg" alt="ar1" width="600" height="300" class="center" />


To generate samples, we define the following function:

{% highlight r linenos %}
ar.samps <- function(n) {
   # empty vector to store results
   res <- numeric(n)
   a.rate <- numeric(n)
   # start loop
   for (t in 1:n) {
      # empty NA object
      tmp.res <- NA
      # counter 
      iter <- 0
      # loop until u <= r
      while (is.na(tmp.res)) {
         # sample h(theta)
         x <- rexp(1, rate=1/2)
         # sample uniform(0,1)
         u <- runif(1)
         # evaluate
         r <- dchisq(x, 2)/(2 * dexp(x, 1/2))
         # accept if u <= r
         if (u <= r) tmp.res <- x
         # update counter
         iter <- iter + 1
      }
      # if accepted store result
      res[t] <- tmp.res
      a.rate[t] <- 1/iter
   }

   message(paste0('Acceptance rate :',round(mean(a.rate),3)))
   return(res)
}
{% endhighlight %}

Now, let's see whether it works.

{% highlight r linenos %}
# set seed
set.seed(123)
# generate 10000 samples
samps <- ar.samps(10000) 
{% endhighlight %}

<pre><code> ## Acceptance rate :0.699
</code></pre>

We might plot the distribution of the samples and compare it to the true density to evaluate whether the method works as expected.

{% highlight r linenos %}
ggplot(data.table(samps),
       aes(x = samps)) +
    geom_histogram(
        bins = 100,
        aes(y = ..density..),
        col = 'white',
        fill = cols[2],
        alpha = .5,
        boundary = 0
    ) +
    geom_line(
        data = data.table(
            x = seq(min(samps), max(samps), .05),
            y = dchisq(seq(min(samps), max(samps), .05), 2)
        ),
        inherit.aes = F,
        aes(x = x, y = y), 
        col =  cols[1]
    ) +
    theme_bw() +
    labs(x = expression(theta), y = expression(f(theta)))
{% endhighlight %}

<img src="{{ site.baseurl }}/assets/img/arsampling2.jpg" alt="ar2" width="500" height="300" class="center"/>

Close! Note that we might improve the acceptance rate by finding a function {% katex %} g^* {% endkatex %} which better approximates the chi-squared distribution from above. A great strength of accept-reject sampling is that we will obtain independent samples from the target distribution, which is not the case in MCMC methods.

A clear shortcoming of the Accept-Reject sampling method, on the other hand, is that it is hard to apply in situations where {% katex %} \theta {% endkatex %} is multidimensional as it will become increasingly difficult to find a function that will dominate {% katex %} p {% endkatex %} for all {% katex %} \theta \in \Theta {% endkatex %} which is, in turn, proportional to a distribution from which we can easily sample. Also, there are a lot of situations where a function of the form {% katex %} g = k \cdot h {% endkatex %} such that {% katex %} \forall \theta \in \Theta, g \ge p {% endkatex %} is difficult to find.
