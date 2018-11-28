---
layout: post
title: Heterogenity and the S-Shaped Diffusion Curve
date: 2017-09-12 01:42
author: baruuum
comments: true
categories: [Diffusion, Formal Models, Quant Stuff, Quotes, Sociology]
---
In a nice article written by Prof. Payton Young, he shows that a diffusion process where different subgroups adopt a new item at different rates will <em>not</em> show an S-shaped diffusion curve as long as the decisions to adopt the new idea are independent of one another. Following his work, this can be succinctly shown as follows. Consider the diffusion process

$$ \frac{dp(t)}{dt} = \lambda(1-p(t)), $$

where $ \lambda $ is the rate of adoption and $ p(t) $ is the proportion of adopters at time $ t $. Given the initial condition $ p(0)=0 $, the solution to this equation is $ p(t) = 1-e^{-\lambda t} $. Clearly, this process will not show an S-shaped curve over time, as $ \dot p(t)=dp(t)/dt $ is strictly concave in $ t $. Now suppose $ \lambda $ differs across subpopulations. Let $ \lambda \sim \mu $, where $ \mu $ is <em>some</em> distribution function <em>not depending on time</em>. This amounts to the assumption that subgroups differ in their rate of adoption, but that the rate will not change over time. Note that the number of adopters within each subgroup will vary across time as there will be fewer individuals who have not already adopted the idea as time proceeds. Let $ p(\lambda, t) $ be the proportion of adopters in subgroup indexed by $ \lambda $ (note that $ \lambda $ is now a random variable). The proportion of adopters in society as a whole at time $ t $ is then

$$ p(t) =\int p(\lambda, t)d\mu(\lambda) = 1- \int e^{-\lambda t}d\mu(\lambda). $$

Differentiating this equation twice with respect to $ t $ shows that $ \ddot p(t) = d^2 p(t)/dt^2 <0 $ for all values of $ t $ <em>irrespective of the distribution</em> $ \mu $. On the other hand, for the trajectory to show an S-shaped curve, it has to be the case that $ \ddot p $ is positive at first and then negative. Thus, heterogeneity in adoption rates across subgroups alone is not sufficient to explain the S-shaped diffusion curve.

On the other hand, even simple dependencies on the proportion of other adopters in the population will result in the desired S-shaped curve. For example, setting $ \lambda(t) = \gamma p(t) $ would result in the differential equation,

$$ \frac{dp(t)}{dt} = \gamma p(t)(1-p(t)), $$

i.e., the well-known logistic equation, which is clearly S-shaped.

Of course, the model above has some strong assumptions. Individuals who have once adopted the new item are, for example, assumed not to switch back. Still I think that this is a very interesting result.
