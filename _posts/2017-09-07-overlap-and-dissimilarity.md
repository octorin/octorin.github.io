---
layout: post
title: Coefficient of Overlapping
date: 2017-09-07 03:14
author: baruuum
comments: true
categories: [Measurement, Overlap, Quant Stuff]
---
An interesting measure that I've encountered in my research on polarization is the **coefficient of overlapping** defined as

$$ OVL(f,g) = \int_X \min\{f(x),g(x)\}dx, $$

where $$ f $$ and $$ g $$ are two density functions defined over the same region $$ X $$. The measure is quite intuitive in that it measures the overlapping region of two densities. If the densities are defined over disjoint subsets of $$ X $$, say $$ X_f $$ and $$ X_g $$, $$ X_f \cap X_g = \varnothing $$, then $$ OVL(f,g)=0 $$. Also, as the area under each of the densities has to integrate to one, $$ OVL(f,g)=1 $$ if and only if the two functions are equal almost everywhere. For the discrete case, we might write

$$ OVL(f,g) = \sum_{x\in X} \min\{f(x),g(x)\}, $$

with $$ f,g $$ being PMFs.

Another interesting measure, encountered in the literature on neighborhood segregation, is the dissimilarity index. Consider two populations $$ A $$ and $$ B $$, which are distributed over the same area. Suppose the area is divided into cells (e.g., census tracks) which we denote by $$ x\in X $$. The dissimilarity index is defined as

$$ D(f,g) = \frac{1}{2}\sum_{x\in X} \left|f(x)-g(x)\right|, $$

where $$ f(x)=n(x)/N $$ and $$ g(x)=m(x)/M $$ with $$ N $$ and $$ M $$ being, respectively, the total number of individuals in group $$ A $$ and $$ B $$, and $$ n(x) $$ and $$ m(x) $$ are the number of individuals from group $$ A $$ and $$ B $$ which are found in cell $$ x $$.

Now, let us go back, for a moment, to the coefficient of overlapping. It is not difficult to see that for two non-negative real numbers $$ a,b\in \mathbb R_+ $$,

$$ \min\{a,b\} = \frac{1}{2}(a+b - |a-b|). $$

So rewriting $$ OVL(f,g) $$ using this formula, we obtain

$$ OVL(f,g) = \int_X \frac{1}{2}\Big(f(x)+g(x)-|f(x)-g(x)|\Big)dx = 1 - \int_X|f(x)-g(x)|dx, $$

with discrete analog

$$ OVL(f,g) = 1-\sum_{x\in X}|f(x)-g(x)|. $$

Thus, it turns out that $$ D(f,g) = 1-OVL(f,g) $$.

Another interesting property of $$ OVL(f,g) $$ is the following. Note that

$$ OVL(f,g) =\int_X \min\{f,g\}dx = \int_X \mathbb I(f<g)dF + \int_X\mathbb I(g<f)dG  $$

Thus, we have

$$ OVL(f,g) = E_F[\mathbb I(f(X)<g(X))] + E_G[\mathbb I(g(X)<f(X))], $$

where $$ E_H(W) $$ is the expectation of the random variable $$ W $$ under the distribution $$ H $$. Thus, the coefficient of overlapping shows the error rate when we would infer the true distribution of an observed data point $$ x $$, amongst the two candidate distributions $$ f $$ and $$ g $$, based on the decision rule of choosing the distribution with a higher density at each point.
