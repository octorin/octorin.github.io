---
layout: post
title: The Random Walk Metropolis Algorithm
date: 2017-12-05 08:29
author: baruuum
comments: true
categories: [Bayesian, MCMC, Metropolis, Monte Carlo Methods, Quant Stuff, Sampling]
---
In a [previous post](/blog/2017/accept-reject-sampling/), we looked at the Accept-Reject sampling method. This method can be used in  situations where we are able to evaluate the distribution function $ p  $ from which we want to sample, but sampling from $ p  $ directly is difficult or inefficient. As noted in that post, the Accept-Reject sampling method is quite easy to implement but has the shortcoming that we have to find a density function $ h  $ form which it is easy to sample and that satisfies the condition that the target density $ p  $ is strictly dominated by the function $ k\cdot h  $, where $ k  $ is a constant. This will become increasingly difficult as the target distribution becomes more complex and the dimensionality of the parameter vector increases. On the other hand, the Metropolis algorithm does not require us to find such a density $ h  $. Instead, it formulates a set of "jumping" rules that form a Markov chain which is guaranteed, under reasonable assumptions, to converge to a unique stationary distribution that conincides with the target density $ p  $ from which we want to sample. Thus, after the chain reaches the stationary distribution, every "step" of the chain generates a sample from the distribution $ p $.

Similarly to Accept-Reject sampling, the basic requirement for the Metropolis algorithm to work is that we can evaluate a function $ f  $ which is proportional to $ p  $. This property is especially important in the Bayesian context, as it is always the normalization constant of the posterior distribution that is most difficult to compute. Yet, the Metropolis algorithm, as all other Markov chain Monte Carlo methods, suffers from the shortcoming that the sampled values will not be independent samples. As it will become clear later, each sampled value depends on the previously sampled value which induces autocorelation between the samples.

<h3>The Metropolis Algorithm</h3>

The Metropolis algorithm can be summarized as follows. Starting from a parameter value $ \theta^t  $, a proposal for the next move is sampled from a <em>symmetric</em> "proposal" or "jumping" distribution $ J(\theta^*,\theta^t)  $. Thereafter a procedure similar to accept-reject sampling is used to decide whether to accept the new value or not. Schematically, the Metropolis algorithm looks something like the following:

<ol>
<li>sample $ \theta^*  $ from $ J(\theta^*,\theta^t)  $.</li>
<li>calculate $ \rho= \min[1, p(\theta^*)/p(\theta^t)].  $</li>
<li>sample $ u \sim \text{unif}(0,1)  $.</li>
<li>if $ u\le \rho$ then $ \theta^{t+1} = \theta^*  $, else set $ \theta^{t+1}=\theta^t  $.</li>
</ol>

Describing why this algorithm belongs to the family of Markov chain Monte Carlo methods will make also clear why it works.

First, given the current value $\theta^t$ and a proposal $\theta^\*$, the algorithm will move to $\theta^\*$ whenever $p(\theta^\*) > p(\theta^t)$. Otherwise, the algorithm moves to $\theta^\*$ with probability $\rho$ and stays at its current value with probability $1-\rho$. Note that the algorithm is ignorant of the history of previous visits and depends solely on where it is right now in determining the next sampled value. Thus, the probabilistic movement of the sampler can be conceptualized as a first-order Markov chain where the <em>transition probabilities</em> are defined in terms of $ J(\cdot,\cdot)  $ and $ \rho$.

For instance, suppose that the distribution from which we want to sample is supported on three finite values, i.e., $\Theta= \\{\theta_1,\theta_2,\theta_3\\}$ with $ p(\theta_i)=p_i, i=1,2,3.  $ Suppose further that the proposal distribution $ J  $ proposes the two states at which we are currently <em>not</em> with equal probability. That is, if we are at $ \theta_2  $, $ J(\theta_1,\theta_2)=J(\theta_3,\theta_2)=1/2  $. Pick any two states $ i  $ and $ j  $ and consider the ratio of moving from $ i  $ to $ j  $ and from $ j  $ to $ i  $, i.e., $ \Pr[\theta_i\rightarrow \theta_j]/\Pr[\theta_j\rightarrow \theta_i]  $. The numerator can be expressed as the probability that $ \theta_j  $ is proposed times the probability that it is accepted given that it was proposed. Thus, we have

$$  \Pr[\theta_i\rightarrow \theta_j] = \frac{\rho}{2} = \frac{\min(p_j/p_i,1)}{2}. $$

$ \Pr[\theta_j\rightarrow\theta_i]  $ can be derived in the same way. Now, it turns out that the ratio of moving to the two different states is proportional to the ratio of the corresponding probabiliites of our target distribution. Namely,

$$ \begin{aligned}\frac{\Pr[\theta_i\rightarrow \theta_j]}{\Pr[\theta_j\rightarrow \theta_i]}&=\frac{.5\times \min(p_j/p_i,1)}{.5\times \min(p_i/p_j,1)}\\
&=\frac{p_j/p_i}{1}\mathbb I(p_i>p_j)+ \frac{1}{p_i/p_j}\mathbb I(p_j>p_i) + 1\cdot \mathbb I(p_i=p_j)\\&=\frac{p_j}{p_i},\end{aligned} $$

where $ \mathbb I(z)  $ is an indicator function that is equal to one if $ z  $ is true and zero otherwise. Thus, the <em>relative</em> probability of moving from a state to another matches exactly the ration of their probability under the target density. This means by sampling the trajectory of the Markov chain, we will be sampling proportional to the target distribution!

Let us make this intuition a little bit more rigorous. We will assume that the parameter space $ \Theta  $ is both discrete and finite. I call it the <em>parameter space</em> as I have a Bayesian model in mind, but this is nothing but the state-space in Markov chains and could be also (as in the example from above) be the support of any random variable. For uncountably infinite parameter spaces, the intuition of what follows will remain the same, although a rigorous treatment would require knowledge of measure theory. For example, if $ \Theta  $ is an uncountable set, the probability to move to a particular $ \theta\in \Theta  $ will be always zero, so that the transition probabilities are ill-defined in the current framework. Instead, we have to conceptualize the transition probabilities in terms of subsets of $ \Theta  $, $ A\in \mathcal F  $, where $ \mathcal F  $ is usually taken as the Borel Sigma-algebra of $ \Theta  $. Here, we will not deal with these complexities but simply assume that all the results on finite parameter spaces carry over to infinite ones.

We assume that for all $ \theta_1,\theta_2\in \Theta\times\Theta,  $ we have $ J(\theta_1,\theta_2) > 0  $. In words, this means that regardless what our current parameter value is, we have a non-zero probability to jump to <em>any</em> other value in the parameter space in the next move. Thus, the condition ensures that we can jump from any state to every other state at each iteration of the Markov chain. This implies, in turn, that the Markov chain constructed by $ J  $ and $ \rho  $ is <em>irreducible</em> (all states in $ \Theta  $ belong to the same communicating class), <em>positive recurrent</em> (for all $ \theta_i\in \Theta  $ the expected waiting time for a chain starting at $ \theta_i  $ to return to state $ \theta_i  $ is finite), and <em>aperiodic</em> (for at least one $ \theta_i\in \Theta  $, we have $ \Pr[\theta_i \rightarrow \theta_i]>0  $. But as all states of $ \Theta  $ belong to the same communicating class, the period of the Markov chain is equal to one). Given these results, standard Markov chain theory tells us that a <em>unique</em> stationary distribution $ \pi  $ of the Markov chain exists. In plain words, this means that the Markov chain will converge to the <em>same</em> distribution $ \pi  $ <em>irrespective of</em> the specific initial point from which we have started the chain.

Our next move is to show that the unique stationary distribution $ \pi  $ of the Markov chain is equal to $ p  $, our target distribution. To sketch the proof, suppose we start the chain at the target distribution $ p  $. The probability of moving from $ \theta^t  $ to $ \theta^{t+1}  $ is

$$  \Pr[\theta^t \rightarrow \theta^{t+1}] = J(\theta^{t+1},\theta^t)\rho = J(\theta^{t+1},\theta^t)\min\left(\frac{p(\theta^{t+1})}{p(\theta^t)},1\right) $$

Note that this defines the <em>transition kernel</em> (or the <em>transition matrix</em> in the finite case) of the Markov chain generated by the algorithm. That is, $ \Pr[\theta_i \rightarrow \theta_j] = (\mathbf {P})_{ij}  $, where $ \mathbf P $ is a row-stochastic matrix. Let us write $ \Pr[\theta^t \rightarrow \theta^{t+1}] = P(\theta^t, \theta^{t+1}).  $ Also, notice that $ p  $ needs only to be proportional to the true target density. For example, if we use $ h(w) = k\cdot p(w)  $, with $ k  $ being a constant, $ h(\theta^{t+1})/h(\theta^t) = p(\theta^{t+1})/p(\theta^t)  $, and, thus, the probability of moving from $ \theta^t  $ to $ \theta^{t+1}  $ is not affected by $ k  $.

Now, multiplying both sides by $ p(\theta^t)  $, we obtain

$$  p(\theta^t)P(\theta^t, \theta^{t+1}) = J(\theta^{t+1},\theta^t)\min\left(p(\theta^{t+1}),p(\theta^t)\right).  $$

But, as the proposal distribution is symmetric, we have also that

$$  p(\theta^{t+1})P(\theta^{t+1},\theta^t) = J(\theta^{t+1},\theta^t)\min\left(p(\theta^{t+1}),p(\theta^t)\right)  $$

so that

$$  p(\theta^t)P(\theta^t, \theta^{t+1}) =p(\theta^{t+1})P(\theta^{t+1},\theta^t).  $$

The last equation is called the <em>detailed balance equation</em> and a Markov chain with stationary distribution satisfying the detailed balance condition is said to be <em>reversible</em>. The important point here is that a distribution $p$ satisfies the detailed balance condition only if $p$ is indeed the stationary distribution of the Markov chain. 

To illustrate, suppose the chain has reached $p$ at its $ t  $th iteration and consider the probability that it will move to state $ \theta_j  $ in the next iteration. This probability can be written as

$$  \sum_{i}p(\theta_i) P(\theta_i, \theta_j)=\sum_{i}p(\theta_j)P(\theta_j,\theta_i)=p(\theta_j)\sum_{i}P(\theta_j,\theta_i) = p(\theta_j).  $$

for all $ j  $, where we used reversibility at the second step and the sum of $P(\theta_j, \theta_i) = \Pr[\theta_j \rightarrow \theta_i]$ over $i$ is equal to one. Thus, if the chain is at $p$, the probability of moving to $\theta_j$ in the next step is equal to $p(\theta_j)$; once the Markov chain has reached this distribution, it will stay there. This is precisely what is meant by a distribution being stationary. Yet, we know that the stationary distribution is unique, which shows that $p$ has to be the distribution we are seeking.

In sum, by irreducibility, positive recurrence, and aperiodicity, the Markov chain constructed by the Metrpolis algorithm has a <em>unique</em> stationary distribution; and further using the detailed balance condition, we have shown that this unique stationary distribution coincides with the target density from which we want to sample. Thus every new value of $ \theta\in\Theta  $ that is visited after the Markov chain has reached the stationary distribution will be a sample from $ p  $.

However, although we know that the chain will <em>eventually</em> converge to $ p  $, it will take some time until the chain reaches the stationary distribution. This is the reason why MCMC methods usually drop several iterations before starting sampling (the so-called "burn-in" period).

<h3>Random-Walk Metropolis in R</h3>

Now, let us try the algorithm out. Consider a situation in which we know in advance that our parameter vector $ \theta=[\theta_1,\theta_2]  $ is distributed according to a bivariate Normal distribution, $ \text{Normal}_2(\mu,\Sigma)  $, but have no way to sample from it. Suppose that $ \mu_1=0,\mu_2=1  $, $ \text{diag}(\Sigma) = \mathbf I_2  $, and the off-diagonal element of $ \Sigma $ is equal to $ 0.5  $. To use the Metropolis algorithm, we have to define a proposal distribution, which is symmetric. Define $ J(\theta^{t+1},\theta^t) = (\theta_1^t + Z_1, \theta_2^t + Z_2)  $ where $ Z_i \sim \text{Normal}(0,\delta^2), i=1,2  $. That is, the proposal distribution is a bivariate Normal distribution with standard deviations $ (\delta,\delta)  $ centered at current parameter value (and zero correlations between the dimensions). 

To calculate $ \rho  $, we have to evaluate some function that is proportional to the target density. For the bivariate normal density with the parameters as defined above, we have

$$ \begin{aligned}
f(\theta|\mu,\Sigma) &=\propto \exp\left(-\frac{2}{3}\Big[\theta_1^2 + (\theta_2-1)^2 - \theta_1(\theta_2-1)\Big]\right)\\
&= h(\theta)
\end{aligned}$$

In coding the algorithm, we use the log of $ h(\theta)$, for numerical stability. First, we load some libraries.

{% highlight r %}
library('data.table')
library('dplyr')
library('ggplot2')
library('ellipse')
library('viridis')
library('gridExtra')
library('ellipse')
{% endhighlight %}

The log of the target density might be coded as follows:

{% highlight r %}
# logged bivariate normal density (up to proportionality constant)
target.fun <- function(theta) {
   -(2/3) * (theta[1]^2 + (theta[2] - 1)^2 - theta[1] * (theta[2] - 1))
}
{% endhighlight %}

Now, we define the Metropolis algorithm.

{% highlight r %}
metropolis <- function(theta,
                       delta = .3,
                       n.iter = 5000,
                       burnin = 2000) {
   # dimensions of theta
   d <- length(theta)
   # object to store results
   res.mat <- matrix(NA, ncol = d, nrow = n.iter - burnin + 1) 
   if (burnin > 0) {

      for (ii in 1:burnin) {

         # J distribution
         u <- rnorm(d, 0, delta)
         # proposal value
         proposal <- theta + u
         # evaluate log-ratio
         lr <- target.fun(proposal) - target.fun(theta)
         if (log(runif(1)) <= lr) theta <- proposal

      }
   }

   res.mat[1,] <- theta

   for (ii in 2:(n.iter - burnin + 1)) {

      # symmetric J distribution
      u <- rnorm(d, 0, delta)
      # proposal value
      proposal <- theta + u
      # evaluate log-ratio
      lr <- target.fun(proposal) - target.fun(theta)
      if (log(runif(1)) <= lr) theta <- proposal
      res.mat[ii,] <- theta

   }

   return(res.mat)

}
{% endhighlight %}

Let us have a look at how the sampler behaves. We define four different starting points $ (\pm 4, \pm 4)  $, and let the sampler run for 50,100, 250, and 500 iterations. We might plot the trajectory of the sampler, together with the region that contains 95% of the density under the target distribution, to see how well it behaves.

{% highlight r %}
# initial values
inits <- list(c(-4, -4), c(-4, 4), c(4, -4),c(4, 4))

# create 95% region of target density
mu = c(0,1)
sigma = matrix(c(1, .5, .5, 1), byrow = T, nrow = 2)
norm.ell <- ellipse(sigma, scale = c(1, 1), centre = mu, level = .95) %>%
            as.data.frame

# how many iterations?
iters <- c(50,100,250,500)

# run sampler
g.list <- vector('list',length(iters))
cc <- 0
for (iter in iters) {
   cc <- cc+1

   # run algorithm
   res <- lapply(1:length(inits), function(w) { 
                  data.table(samples=metropolis(theta = inits[[w]], 
                                                delta = .5, 
                                                burnin = 0, 
                                                n.iter = iter), 
                             chain.no = w) %>%
                  setnames(c('theta.1', 'theta.2', 'chain'))
               }
         )

   # combine results into data.table
   df <- rbindlist(res)

   # save plot
   g.list[[cc]] <- ggplot(df, 
                      aes(x = theta.1,y = theta.2, group = chain)) +
                   geom_path(col = viridis(1), alpha = .5) +
                   geom_point(data = as.data.frame(
                                 cbind(do.call('rbind', inits),
                                       1:4)),
                              aes(x = V1, y = V2, group = V3)) +
                   geom_path(data = norm.ell, inherit.aes = F,
                             aes(x = x, y = y),
                             size = .7)+
                   theme_bw()+
                   scale_x_continuous(breaks = seq(-6, 6, 2)) +
                   scale_y_continuous(breaks = seq(-6, 6, 2)) +
                   labs(x = expression(theta[1]), 
                        y = expression(theta[2]))+
                   coord_fixed() +
                   theme(plot.margin = grid::unit(rep(0, 4), 'mm'))+
                   ggtitle(paste0('Iterations = ', iter))
}

# plot results
grid.arrange(grobs = g.list, nrow = 2)
{% endhighlight %}


<center><img src="/assets/img/metropolis1.jpg" width="500" height="500" /></center>

As it is clear from the plot, the sampler needs some time to find the stationary distribution. But once it has converged, it will sample from the desired target density. Thus, to retain only the "valid" samples, we will have to discard those sampled values at the start of the iterations.

Another important point to note is that the performance of the sampler will crucially depend on the propsal distribution, $ J  $. Recall that $ \delta  $ is the parameter that tunes the variance of the proposal distribution in our algorithm, i.e., $ \theta_i^* \sim \text{Normal} (\theta_i^t,\delta^2), \quad i=1,2$, where $ \theta_i^*$ is the proposed next step and $ \theta_i^t  $ is the current value that was sampled. If $ \delta  $ is too small, the proposed steps will be concentrated around where the sampler finds itself. Therefore, it will have a hard time to find the stationary distribution. And even if the stationary distribution is found, the sampler is inefficient in exploring the whole parameter space as the next step will very close to the current one. On the other hand, if $ \delta$ is too large, the proposed values will be spread around the whole space. This implies that proposals at low-density region of the target distribution will be proposed quite often, so that too many proposals will be rejected, which leads again to an inefficient behavior. To illustrate, let us run the sampler again, but this time fixing the number iterations to 500 and changing the value of $ \delta  $.

{% highlight r %}
# different values for delta
d.vec <- c(.05, .1, 1, 3)
# run sampler
g.list <- vector('list', length(iters))
cc <- 0
for (d in d.vec) {
   cc <- cc + 1

   # run algorithm
   res <- lapply(1:length(inits), function(w) {
                  data.table(samples = metropolis(theta = inits[[w]], 
                             delta = d, burnin = 0, n.iter = 500), 
                             chain.no = w) %>%
                  setnames(c('theta.1', 'theta.2', 'chain'))
                 }
          )
   # combine results into data.table
   df <- rbindlist(res)

   # save plot
   g.list[[cc]] <- ggplot(df, 
                     aes(x = theta.1, y = theta.2, group = chain)) +
                     geom_path(col = viridis(1), alpha = .5) +
                     geom_point(data = as.data.frame(cbind(
                                do.call('rbind', inits), 1:4)),
                                aes(x = V1, y = V2, group = V3)) +
                     geom_path(data = norm.ell, inherit.aes = F,
                               aes(x = x, y = y),
                               size = .7)+
                     theme_bw()+
                     scale_x_continuous(breaks = seq(-6, 6, 2)) +
                     scale_y_continuous(breaks = seq(-6, 6, 2)) +
                     labs(x = expression(theta[1]), 
                          y = expression(theta[2]))+
                     coord_fixed() +
                     theme(plot.margin=grid::unit(rep(0, 4), 'mm'))+
                     ggtitle(paste0('Delta = ', d))
}

# plot results
grid.arrange(grobs = g.list, nrow = 2)
{% endhighlight %}

<center><img  src="/assets/img/metropolis2.jpg" width="500" height="500" /></center>

As expected, for too small values of $ \delta  $, the sampler is very inefficient, while for too large values, the sampler takes much fewer steps; and if it moves, it moves in large steps. For the Normal distribution, analytical results exist that determine the optimal values of $\delta$. Yet, for other distributions, the only way to choose a good $\delta$ is often to monitor the acceptance rates, i.e., how often a proposal was accepted/rejected. If the acceptance rate is too high, this means that the sampler moves in too small steps, so increasing $\delta$ would give better results; the opposite holds for the case in which the acceptance rate is too low.

A last important point, which was also noted above, is that the Metropolis algorithm as other MCMC mdethods will generate generate <em>dependent</em> samples. Indeed the next sampled value will inherently depend on where the sampler was at the previous step, which is clear from the <em>traceplot</em> (shown below). But even from the figure above, especially for small values of $ \delta  $, we can see immediately that the samples are not independent. The autocorrelation between subsequent samples can be reduced by sampling only the $k$th step of the sampler, which is called <em> thinning</em>. Yet, using more efficient algorithms, such as Hamiltonian Monte Carlo, is probably a better way to go, especially if all the parameters in the model are continuous.

{% highlight r %}
res <- metropolis(theta = c(-4, 4), 
                       delta = 1,
                       burnin = 0,
                       n.iter = 2000)
df <- melt(res) %>%
      setnames(c('iter', 'var', 'val')) %>%
      data.table

ggplot(df, aes(x = iter, y = val)) + 
   geom_line(col = viridis(1), alpha = .8) + 
   theme_bw() +
   facet_wrap(~var, nrow=2) +
   labs(y = 'Sampled Values',x = 'Iteration') +
   ggtitle(expression('Traceplot of'~theta))
{% endhighlight %}

<center><img  src="/assets/img/metropolis3.jpg" width="500" height="400" /></center>
