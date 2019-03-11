---
layout: post
title: The Spatial Voting Model
date: 2017-08-29 02:53
author: baruuum
comments: true
categories: [Formal Theory, Political Science, Quant Stuff]
---
One interesting model in the political science literature is the spatial voting model. The basic idea behind it is to assume that each voter (or legislator) has an "ideal point" in a (possibly multidimensional) policy space and votes for the policy proposal closest to this point. It turns out that if we assume that voters are have quadratic utility functions with Normally distributed stochastic errors, the resulting model is closely related to the two-parameter Normal ogive IRT model in the psychology literature.

To show how this model works, suppose that on the $$ j $$th policy proposal, voter $$ i $$ chooses between two alternative $$ \psi_j \in \mathbb R^d $$ and $$ \zeta_j \in \mathbb R^d $$, where $$ \psi_j $$ is the position of the "yea" outcome and $$ \zeta_j $$ is the "nay" outcome. The utility that voter $$ i $$ receives from each of the alternatives is assumed to be

$$ 
u_i(\psi_j) = -\frac{1}{2}\|\theta_i-\psi_j\|^2 + \epsilon_{ij}^\psi  \qquad \text{and} \qquad u_i(\zeta_j) = -\frac{1}{2}\vert\theta_i-\zeta_j\vert^2 + \epsilon_{ij}^\zeta 
$$

where $$ \theta_i \in \mathbb R^d $$ is the voter's "ideal point," $$ \vert\cdot\vert $$ is the Euclidean norm, and the $$ \epsilon $$'s are stochastic errors. The voter prefers the outcome $$ \psi_j $$ over $$ \zeta_j $$ if

$$ u_i(\psi_j) > u_i(\zeta_j) $$

which implies that

$$ 
\epsilon_{ij}^\zeta - \epsilon_{ij}^\psi < \frac{1}{2}(\vert\zeta_j\vert^2 - \vert\psi_j\vert^2) + (\psi_j-\zeta_j)^\top \theta_i = \alpha_j + \beta_j^\top \theta_i,
$$

where $$ \alpha_j \in \mathbb R $$ and $$ \beta_j \in \mathbb R^d $$. If $$ \epsilon_{ij}^\zeta - \epsilon_{ij}^\psi \sim G $$, where $$ G(x) = \text{Prob}[\epsilon_{ij}^\zeta - \epsilon_{ij}^\psi \le x], x\in \mathbb R $$, it follows that the probability that the voter prefers the outcome $$ \psi_j $$ over $$ \zeta_j $$ is $$ G(\alpha_j + \beta_j^\top \theta_i) $$.

Thus, if we let $$ G $$ be a Normal distribution, then we get the two-parameter normal ogive model; if $$ G $$ is assumed to be a Logistic distribution, we have the usual two-parameter IRT model. 
