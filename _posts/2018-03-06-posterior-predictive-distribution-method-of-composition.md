---
layout: post
title: Posterior Predictive Distribution & Method of Composition
date: 2018-03-06 08:53
author: baruuum
comments: true
categories: [Bayesian, Method of Composition, Monte Carlo Methods, Posterior Predictive Distribution, Quant Stuff]
---

**DISCLAIMER**: _This is a blog post from my graduate school years. The blog is no longer maintained and the material might include typos and errors._


The posterior predictive distribution (PPD) of unobserved (or future) data are of major importance in Bayesian model checking or imputation methods that rely on a Bayesian framework. At first sight, the connection between computationally obtaining samples from the PDD and the mathematical expression of it are not obvious. So here is why it works.

Let $ y $ be a vector of observed data that come from a probability distribution $ p(y\vert\theta) $ and let $ \tilde y $ be a vector of future (or out-of-sample) values we want to predict. We assume that $ \tilde y $ comes from the same distribution as $ y $, so the distribution of interest is $ p(\tilde y \vert y) $, i.e., the distribution of out-of-sample values given our observed sample. It might be tempting to use our best estimate---such as the MLE or MAP estimate---to obtain information about the PPD of future data. However, doing so would inevitably ignore our uncertainty about $ \theta $. Thus, the appropriate way to procede is to average over the posterior distribution of $ \theta $, namely $ p(\theta\vert y) $. Notice also that $ \tilde y $ is independent of $ y $ given $ \theta $, as it is assumed to be an independent sample drawn from the same distribution as $ y $. Thus,

$$
p(\tilde y\vert \theta, y) = \frac{p(\tilde y, y\vert\theta )p(\theta)}{p(\theta, y)} =
\frac{p(\tilde y\vert\theta )p(y \vert\theta) p(\theta)}{p(y\vert \theta)p(\theta)} = p(\tilde y \vert\theta).
$$

The posterior predictive distribution of $ \tilde y $ is thus,

$$  p(\tilde y\vert y ) = \int_\Theta p(\tilde y \vert \theta,y) p(\theta \vert y) d\theta = \int_\Theta p(\tilde y \vert \theta) p(\theta \vert y) d\theta  $$

where $ \Theta $ is the support of $ \theta $.

Now, how do we obtain the distribution $ p(\tilde y\vert y) $? While the evaluation of the integral seems to be complex, we can easily appximate it by drawing samples from $ p(\tilde y \vert y) $ given that we have a set of posterior draws from $ p(\theta\vert y) $. That is,

**for** $s$ in $1:S$ **do**                                   
 1. draw $ \theta^{(s)} $ from $ p(\theta\vert y) $               
 2. draw $ \tilde y^{(s)} $ from $ p(\tilde y\vert\theta^{(s)}) $ 
 
where, in most situations, we have already the draws from $ p(\theta\vert y) $, so that only the second step is required.

The reason why this works is based on a simulation method known as the <em> method of composition </em>.

Considering that $ p(\tilde y, \theta \vert y) = p(\tilde y\vert \theta, y)p(\theta \vert y) $, first sampling a parameter vector $ \theta^{(s)} $ from the posterior, $ p(\theta\vert y) $, and then using this vector to sample $ \tilde y^{(s)} $ from the perdictive distribution, $ p(\tilde y \vert \theta) = p(\tilde y \vert \theta, y) $, gives us a sample $ (\tilde y^{(s)}, \theta^{(s)}) $ from the joint distribution $ p(\tilde y, \theta\vert y) $. It follows that the sampled values $ \tilde y^{(s)}, s=1,2,...,S $, are samples from the posterior predictive distribution $ p(\tilde y \vert y) $.
