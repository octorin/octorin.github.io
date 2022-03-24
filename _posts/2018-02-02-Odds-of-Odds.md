---
layout: post
title: Ratio of Odds Ratios
date: 2018-02-02 08:23
author: baruuum
comments: true
categories: [Logistic Regression, Odds Ratios, Quant Stuff]
---

**DISCLAIMER**: _This is a blog post from my graduate school years. The blog is no longer maintained and the material might include typos and errors._


Logistic regression is probably one of the most often used models in sociology. The general logistic regression model might be written as

$$  \ln\left(\frac{p(\mathbf{x})}{1-p(\mathbf{x})}\right) = \mathbf{x}^\top \boldsymbol{\beta} $$

 where

 $$  p(\mathbf{x}) = E[y \vert\mathbf{x}]  $$

and $ y $ is a binary outcome, $ \mathbf{x} $ a $ k $-dimensional vector of predictors, and $ \boldsymbol{\beta} $ a vector of coefficients. In the last class in which I am TAing, a student asked about the interpretation of exponentiated logistic regression coefficients in terms of odds ratios.

Of course, my first response was to say that they should always plot the predicted probabilities from their models, since nobody really understand what odds-ratios are. But then, I became myself curious what the exponentiated interaction coefficient in a logistic regression represents.

So, suppose we have two dummy variables and their interaction in our model:

$$  \ln\left(\frac{p(\mathbf{x})}{1-p(\mathbf{x})}\right) = \beta_0 + \beta_1x_1 + \beta_2x_2 + \beta_3x_1x_2. $$

By exponentiating both sides we obtain

$$  Odds(x_1,x_2) = \frac{p(\mathbf{x})}{1-p(\mathbf{x})} = \exp(\beta_0 + \beta_1x_1 + \beta_2x_2 + \beta_3x_1x_2). $$

By pluggin-in different values of $ x_1 $ and $ x_2 $, we therefore have

$ Odds(x_1=0,x_2=0) = \exp(\beta_0) = \gamma_0 $
$ Odds(x_1=1,x_2=0) = \exp(\beta_0 + \beta_1) = \gamma_0\gamma_1 $
$ Odds(x_1=0,x_2=1) = \exp(\beta_0 + \beta_2) = \gamma_0\gamma_2 $
$ Odds(x_1=1,x_2=1) = \exp(\beta_0 + \beta_1 + \beta_2 + \beta_3) = \gamma_0\gamma_1\gamma_2\gamma_3 $

Hence, it follows that

$$
\begin{aligned}
\exp(\beta_3) &= \gamma_3 = \frac{\gamma_0\gamma_1\gamma_2\gamma_3\gamma_0}{\gamma_0\gamma_1\gamma_0\gamma_2} = \frac{Odds(x_1=1,x_2=1)Odds(x_1=0,x_2=0)}{Odds(x_1=1,x_2=0)Odds(x_1=0,x_2=1)} \\
&=\frac{Odds(x_1=1,x_2=1)/Odds(x_1=0,x_2=1)}{Odds(x_1=1,x_2=0)/Odds(x_1=0,x_2=0),}
\end{aligned} 
$$

which is, oddly enough, a ratio of odds-ratios....
