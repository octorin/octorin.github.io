---
title : Regualrity Conditions and MLE
date : 2020-04-25 00:00:00
author : baruuum
comments : true
---

**DISCLAIMER**: _This is a blog post from my graduate school years. The blog is no longer maintained and the material might include typos._

In my first year of graduate school, I set a couple of goals that I wanted to achieve before graduating. One of them was understanding Harrison C. White's *Identity and Control*, which felt, as Randall Collins put it somewhere, like an "IQ test for sociologists'' (6 years later, I am still not sure whether I understood that book). Another goal that I set for myself was understanding the method of maximum likelihood estimation (MLE). Using compiled packages to estimate GLMs via MLE or coding-up Newton-type algorithms was fairly easy, as both of these tasks are quite mechanical. But the method per se remained a mystery to me for quite long, mainly because I was not able to understand the **regularity conditions** that are always mentioned but not often explained. So, sometime ago, I started to write a document that summarizes the regularity conditions in my own words. Much of the details come from Lehmann and Casella (1998)'s [book on the theory of point estimation](https://link.springer.com/book/10.1007/b98854).

As I have decided to close this blog in the near future, I thought this topic would be great for the last post since my graduate years are coming to an end as well. The document basically tries to explain why all these conditions are needed and points out where they appear in the proves of major theorems related to MLE. As always, equations are unnecessarily lengthy, not to clutter the document, but to make them more approachable (i.e., so that one can just "read'' through them). In the end, however, understanding these conditions didn't really influence how I run models in anyway; it was just a small curiosity sparked in my first year of graduate school that kept pushing me deeper into statistical theory. Anyways, here's the document. Most theorems that are mentioned in the document, except the most elementary ones, are explained and proven in the appendix. The only exception is the Radon-Nikodym Theorem.



<h3> Regularity Conditions, Once and For All </h3>



We consider a family of probability distributions $\mathcal P =\\{P\_\theta : \theta\in \Theta\\}$, where $\Theta$ is the parameter space. The parameter space is nothing but the set of values that the parameter $\theta$ can take on, and by a family of distributions we mean a set of distributions that are indexed by the same parameter. For example, we can think of the set of Normal distributions as a family of distributions indexed by $\theta = (\mu, \sigma)$. Each combination of $\mu$ and $\sigma$ will lead to a different distribution, but they would all have the same structure. In the case of the Normal distribution, the parameter space would be $\Theta = \mathbb R \times \mathbb R\_+$ as $\mu$ can be any real number and $\sigma$ any positive real number.

We will assume that the distributions in $\mathcal P$ meet the following conditions:

**Assumption 1**
- **A1** *The distributions $P\_\theta$, $\theta\in \Theta$, are absolutely continuous with respect to a sigma-finite measure $\mu$.*
- **A2** *The parameter space is a subset of $\mathbb R$, i.e., $\Theta \subseteq \mathbb R$.*

I remember that it was intimidating to find out that the very first assumption started with some measure-theoretic technicalities. But the essence of **A1** is quite straightforward: the requirement that the distribution is absolutely continuous with respect to a sigma-finite measure $\mu$ is, by the Radon-Nikodym Theorem, a necessary and sufficient condition for a density $f\_\theta$ corresponding to $P\_\theta$ to exist. So, we can understand **A1** as saying that for each $\theta \in \Theta$, the density $f\_\theta$ exists. Also, in what follows, we will simply write $\int f\_\theta(x) dx$ instead of $\int f\_\theta \mu$ or $\int f\_\theta(x) d\mu(x),$ while keeping in mind that the density $f\_\theta$ is defined always with respect to some base or dominating measure $\mu$. Compared to **A1**, assumption **A2** is very simple and just says that we consider only one-dimensional real-valued parameters. The case in which $\Theta$ is multidimensional, as in the case of the Normal distribution, will require similar but slightly different sets of regularity conditions not discussed in this document.

Next, we define the likelihood function:

**Definition 1** *For a (observed) sample point $\mathbf x = (x\_1,...,x\_n)$ from density $f\_\theta$, the **likelihood function** $L(\theta; \mathbf x) = f\_\theta(\mathbf x)$ is the sample density considered as a function of $\theta$ for fixed $\mathbf x$. The **log-likelihood function** is accordingly defined as $\ell(\theta;\mathbf x) = \log L(\theta;\mathbf x)$.*

We will switch back and forth between the notation $f\_\theta(x)$ and $f(x;\theta)$ when denoting a density with parameter $\theta$. Also, from now on, $\theta\_0$ will always denote the "true'' parameter that has generated the data. The first set of regularity conditions for maximum likelihood estimation are the following:

**Assumption 2 (Regularity Conditions)**
- **R0** *The distributions $P\_\theta$ of the observations are distinct*
- **R1** *The distributions $P\_\theta$ have common support*
- **R2** *The observations are $\mathbf X\_n = \\{X\_1,...,X\_n\\}$, where $X\_i$ are iid with probability density $f\_\theta$ with respect to $\mu$.*

Condition **R0** means that $\theta\ne \theta'$ implies $P\_\theta \ne P\_{\theta'}.$ This condition is necessary to identify the true parameter that generated the data. Suppose, for example, that two different parameters, $\theta\_1$ and $\theta\_2$, lead to the exact same distribution $P$. Then, even if we figure out that the data came from the distribution $P$, we would not be able to tell which of $\theta\_1$ and $\theta\_2$ is the "true'' parameter value. **R1** states that the region on which the distributions are supported does not depend on $\theta$. For continuous densities, this could be understood as saying that the set $\\{x: f\_\theta(x) > 0\\}$ is the same for all $\theta\in \Theta$. Hence, we rule out the case in which, for example, an event $\\{X\_i \le x\_i\\}$ can occur with positive probability when $\theta = \theta\_1$ but not when $\theta = \theta\_2$---i.e., regardless of the value of $\theta$, the random variable $\mathbf X\_n$ will take on the same set of values with positive probability (although what these positive probabilities are will differ according to the value of $\theta$). **R3** says that $\mathbf X\_n$ consists of random variables that are independent and come from the same distribution (the "with respect to $\mu$" is there because all densities are always defined with respect to a dominating measure). The reason why this is important will become clear in the proof of our first theorem.

**Theorem 1** *Under conditions **R0** to **R2**,*

$$P_{\theta_0}[L(\theta_0; \mathbf X_n) > L(\theta; \mathbf X_n)] \longrightarrow 1\text{ as } n\rightarrow\infty$$

*for any fixed $\theta\ne\theta_0$.*

---
We need **R2** to factor the likelihood as

$$L(\theta; \mathbf x_n) = f(\mathbf x_n; \theta) = \prod_{i=1}^n f(x_i;\theta).$$

From this, it follows that

$$
\begin{aligned}
L(\theta_0; \mathbf x_n) > L(\theta; \mathbf x_n) &\iff  \log L(\theta_0;\mathbf x_n) > \log L(\theta;\mathbf x_n) \\
&\iff \sum_{i=1}^n \Big[\log f(x_i; \theta) - \log f(x_i;\theta_0)\Big] < 0  \\
&\iff \frac{1}{n}\sum_{i=1}^n \log\left[\frac{f(x_i;\theta)}{f(x_i; \theta_0)}\right] < 0.
\end{aligned}
$$

By **R1** the ratio $f(x;\theta)/f(x;\theta_0)$ is well-defined for all $x \in \mathcal X$, where $\mathcal X$ is the region of common support. Hence, we can write

$$
\begin{aligned}
\text{E}_{\theta_0}\left[\frac{f(X_i;\theta)}{f(X_i;\theta_0)}\right] &= \int_{\mathcal X} \left[\frac{f(x_i;\theta)}{f(x_i;\theta_0)}\right] f(x_i;\theta_0) dx_i \\
&=\int_{\mathcal X} f(x_i;\theta)dx_i\\ &= 1.
\end{aligned}
$$

Now, by Jensen's Inequality, we have

$$\text{E}_{\theta_0}\log\left[\frac{f(X_i;\theta)}{f(X_i; \theta_0)}\right] < \log\text{E}_{\theta_0}\left[\frac{f(X_i;\theta)}{f(X_i; \theta_0)}\right] =0$$

as the logarithm is a strictly concave function on $\mathbb R_+$. As the expectation on the left-hand side of the inequality is finite, the Weak Law of Large Numbers applies and we have

$$\frac{1}{n}\sum_{i=1}^n \log\left[\frac{f(X_i;\theta)}{f(X_i; \theta_0)}\right] \longrightarrow \text{E}_{\theta_0}\log\left[\frac{f(X_i;\theta)}{f(X_i; \theta_0)}\right],$$

in probability as $n\rightarrow\infty$. But the expectation on the right is a constant smaller than zero. Therefore,

$$\lim_{n\rightarrow\infty} P_{\theta_0}[L(\theta_0; \mathbf X_n) > L(\theta; \mathbf X_n)]  = \lim_{n\rightarrow\infty} P_{\theta_0}\left[\frac{1}{n}\sum_{i=1}^n \log\left(\frac{f(X_i;\theta)}{f(X_i; \theta_0)} \right)< 0\right] =  1,$$

which proves the theorem.
<div style="text-align: right"> Q.E.D </div>
---

This theorem establishes that the likelihood will be strictly larger at the true parameter value than any other value in the parameter space with probability tending to one. Hence, choosing the value $\theta$ that maximizes the likelihood function for any observed dataset seems to be a reasonable choice for a good estimator.

**Definition 2** *Let $\hat\theta(\mathbf x\_n)$ be the value that maximizes the likelihood at the observed $\mathbf X\_n = \mathbf x\_n$. If $\hat\theta(\mathbf X\_n)$ is unique, we call it the **maximum likelihood estimator (MLE)** of $\theta\_0$.*

While Theorem 1 establishes that the likelihood evaluated at $\theta\_0$ will exceed that of any other $\theta\in\Theta$ in the limit, $\hat\theta(\mathbf X\_n)$ might not be unique or might not even exist for finite $n$. In these cases, the MLE is undefined. It turns out that there is a unique condition under which assumptions **R0** to **R2** are sufficient to ensure that the maximum likelihood estimator will be consistent. This is summarized in the following theorem.

**Theorem 2** *If conditions **R0** to **R2** are met and if, in addition, the parameter space is finite, then the MLE exists, is unique with probability tending to one, and is consistent.*

---
As $\Theta$ is finite, we may write $\Theta = \\{\theta_0, \theta_1,...,\theta_K\\}$, where $\theta_0$ is the true parameter. Let

$$\hat\theta(\mathbf X_n) = \text{argmax}_{\theta\in \Theta}L(\theta; \mathbf X_n) = \text{argmax}_{\theta\in\Theta}\ell(\theta;\mathbf X_n)$$

be the MLE, where $\ell$ is the log-likelihood, i.e., $\ell (\theta; \mathbf X_n) = \log L(\theta; \mathbf X_n)$. We want to show that there is a unique value $\hat\theta(\mathbf X_n)$ that maximizes the likelihood and that it tends to $\theta_0$ in probability as $n \rightarrow\infty$.

For each $1\le j\le K$, define the event $A_{jn} =\\{\mathbf X_n : L(\theta_0; \mathbf X_n)>L(\theta_j ; \mathbf X_n)\\}$. Then, Theorem 1 shows that $P_{\theta_0}[A_{jn}] \rightarrow 1$ as $n\rightarrow\infty$ for all $j$. What remains to be shown is that $P_{\theta_0}[A_{1n}\cap A_{2n}\cap \cdots \cap A_{Kn}] \rightarrow 1$ as $n\rightarrow\infty$ as well. It will suffice to show that this result holds for any $A_{jn}\cap A_{j'n}$, $j\ne j'$ (as we can repeatedly apply this result $K-1$ times).

We note that

$$P_{\theta_0}[A_{jn}\cap A_{j'n}] = 1- P_{\theta_0}[A_{jn}^C\cup A_{j'n}^C] \ge 1 - P_{\theta_0}[A_{jn}^C] - P_{\theta_0}[A_{j'n}^C]$$

for all $n$ by the sub-additivity of probability measure. It follows that,

$$ P_{\theta_0}[A_{jn}\cap A_{j'n}] = 1$$

as $n\rightarrow\infty$, since $P_{\theta_0}[A_{jn}^C]$ and $P_{\theta_0}[A_{j'n}^C]$ converge to zero by Theorem 1. Hence,

$$\lim_{n\rightarrow \infty} P_{\theta_0}\left[\bigcap_{k=1}^KA_{kn}\right] = 1.$$

In other words, $\theta_0$ is the unique value at which the likelihood function is simultaneously strictly larger than any other possible parameter value with probability converging to $1$ as $n\rightarrow\infty$. Thus, by choosing the value of $\theta$ that maximizes the likelihood function---i.e., $\hat\theta(\mathbf X_n)$---we will choose the true parameter value $\theta_0$ with probability converging to one. This is equivalent to the statement that

$$P_{\theta_0}[\hat\theta(\mathbf X_n) = \theta_0] \rightarrow 1$$

as $n\rightarrow\infty$. Thus, $\hat\theta(\mathbf X_n)$ is consistent.
<div style="text-align: right"> Q.E.D </div>
---

It turns out that we need stronger assumptions to ensure consistency of the MLE if $\Theta$ is at least countably infinite. While Lehman and Casella show counter-example in which consistency breaks down when $\Theta$ is countable, I couldn't find more general statement (of course, one counter-example suffices to make the case). So, we'll take this as a fact and proceed with infinite parameter spaces.

We add two more regularity conditions

**Assumption 3 (Regularity Conditions, cont)**
- **R3** *The parameter space $\Theta$ contains an open set $\mathcal O$ of which the true parameter value $\theta\_0$ is an interior point.*
- **R4** *For almost all $x\_i$, $f(x\_i;\theta)$ is differentiable with respect to $\theta \in \mathcal O$, with derivative $f'(x\_i;\theta)$.*

**R3** ensures that there is a neighborhood around the true parameter on which the derivative can be defined. **R4** then asserts the existence of the derivative of the density function with respect to $\theta$. The statement "for almost all $x\_i$'' means that a condition holds on a set $A$ that satisfies the condition $P[A] = 1$. In other words, the derivative might not exist for points $x$ that have zero probability of occurrence.

Throughout, we will denote derivatives of the log-likelihood with respect to the parameter $\theta$ as $\ell', \ell\'\',...$. Under the additional assumptions, we have the following result.

**Theorem 3** *Assume conditions **R0** to **R4** hold. Then, with probability tending to 1 as $n\rightarrow\infty$, the* **likelihood equation**

$$\ell'(\theta; \mathbf X_n) = 0$$

*has a root $\hat\theta_n = \hat\theta(\mathbf X_n)$ such that $\hat\theta(\mathbf X_n)$ tends to the true value $\theta_0$ in probability.*

---
By **R3**, there exists $\epsilon>0$ such that $(\theta_0 - \epsilon, \theta_0 + \epsilon) \subset \mathcal O$. Let

$$A_n = \{\mathbf X_n : \ell(\theta_0; \mathbf X_n) > \ell(\theta_0-\epsilon; \mathbf X_n) \text{ and } \ell(\theta_0; \mathbf X_n) > \ell(\theta_0 + \epsilon ; \mathbf X_n)\}.$$

Then, by Theorem 1, $P_{\theta_0}[A_n] \rightarrow 1$ as $n\rightarrow\infty$. Hence, for all $\mathbf X_n \in A_n$, there exists a value $\theta_0 - \epsilon < \hat\theta_n <\theta_0 + \epsilon $ at which $\ell(\theta)$ attains a local maximum, so that $\ell'(\hat\theta_n;\mathbf X_n) = 0$. (The derivative of the log-likelihood equation exists by **R4**. Indeed, it has not even to exist everywhere on $\mathcal X$, but only on a set with probability converging to $1$). This implies that

$$A_n \subseteq B_n, \text{ where }B_n = \{\mathbf X_n : \vert\hat\theta_n(\mathbf X_n) - \theta_0\vert < \epsilon\}\cap \{ \ell'(\hat\theta_n; \mathbf X_n) = 0\}.$$

But for all $n$, $P_{\theta_0}[A_n] \le P_{\theta_0}[B_n]$, and, thus

$$ 1  = \lim_{n\rightarrow\infty}P_{\theta_0}(A_n) \le \limsup_{n\rightarrow\infty}P_{\theta_0}[B_n]\le 1.$$

Hence, for a sequence of roots, $\hat\theta_n$, satisfying the likelihood equations, we have

$$ \lim_{n\rightarrow\infty} P_{\theta_0}\Big[\{\mathbf X_n: \vert\hat\theta_n(\mathbf X_n) - \theta_0\vert<\epsilon\}\Big] = 1$$

for arbitrary small $\epsilon$.

The only problem that remains is that the sequence of roots depend on $\epsilon$. To determine a sequence of roots that does not depend on $\epsilon$, let $\hat\theta_n^{\*}$ be the root "closest'' to $\theta_0$. To see that such a sequence exists, notice that for each $n$, we have $\hat\theta_n \in (\theta_0 - \epsilon, \theta_0 + \epsilon)$. So, the set of roots is a bounded set of real numbers and has an infimum. Therefore, we can choose a sequence of roots $\hat\theta_n^{\*}$, such that for each $n$, $\hat\theta_n^{\*}$ satisfies $\vert\hat\theta_n^{\*} - \theta_0\vert = \inf \\{\hat\theta_n: \vert\hat\theta_n - \theta_0\vert\\}$. If there are two roots that meet this condition, we simply choose the smaller one. This gives an sequence of roots $\hat\theta_n^{\*}$, independent of $\epsilon$, that satisfies

$$\lim_{n\rightarrow\infty}P_{\theta_0}[\vert\hat\theta_n^{*} - \theta_0\vert<\epsilon] \rightarrow 1$$

and we are done.
<div style="text-align: right"> Q.E.D </div>
---

Notice that this theorem establishes only the existence of a sequence of roots to the likelihood equations that converge to the true parameter; it doesn't tell us which sequence of roots to choose as $\theta\_0$ is unknown (we don't know what is the "closest'' root, only that it exists). Of course, if the likelihood equations have a unique root, the MLE will be consistent. Hence, we have the following corollary.

**Corollary 1** *Assume conditions **R0** to **R4** hold. If, in addition, for all $n$ and all $\mathbf x\_n \in \mathcal X$ the likelihood equation has a unique root, $\hat\theta\_n^{\*}$, then the MLE,  $\hat\theta\_n^{\*}$, is a consistent estimator of $\theta\_0$.*

There is also a weaker counter-part:

**Corollary 2** *Assume conditions **R0** to **R4** hold. If, in addition, the probability of the likelihood equation having multiple roots tends to zero as $n\rightarrow\infty$, then the MLE , $\hat\theta\_n^{\*}$, exists for sufficiently large $n$ and is a consistent estimator of $\theta\_0$.*



While Theorem 3 shows that a consistent sequence of roots exist for uncountable parameter spaces, deriving its distribution requires additional regularity conditions.



**Assumption (Regularity Conditions, cont.)**
- **R5** For almost all $x\_i$, $f(x\_i;\theta)$ is twice differentiable with respect to $\theta\in \mathcal O$, and the second derivative is continuous in $\theta$.
- **R6** The integral $\int f(x\_i;\theta)d x\_i$ can be differentiated twice under the integral sign with respect to $\theta\in \mathcal O$.



Notice that **R5** implies **R4**. Condition **R6** is stating that we can freely interchange integration and differentiation in an open neighborhood around $\theta\_0$, i.e.,

$$\frac{\partial^2}{\partial\theta^2}\int_{\mathcal X} f(x_i;\theta) dx_i = \int_{\mathcal X} \frac{\partial^2}{\partial\theta^2}f(x_i; \theta)dx_i.$$

In general, this is possible if the following conditions are met

1. $\text{E}_\theta[(\partial/\partial\theta)f(x\_i;\theta)] < \infty$ for all $\theta\in\mathcal O$.
2. For almost all $x\_i$, the second derivative $(\partial^2/\partial\theta^2)f(x\_i;\theta)$ exists for all $\theta\in\mathcal O$
3. There exists a function $M:\mathbb R \rightarrow \mathbb R$ such that $\vert(\partial^2/\partial\theta^2)f(x\_i;\theta)\vert < M(x)$ for all $\theta\in\mathcal O$ and almost all $x\_i$, where $\text{E}_\theta[M(X\_i)] < \infty$.

Out of these **R5** suffices to ensure Conditions 1. and 2., while Condition 3. states that the second partial derivative is bounded by an integrable function (so we can use the Dominated Convergence Theorem).

**Lemma 1** *Under condition **R0** to **R6**, the following results hold:*

$$\text{E}_{\theta}\left[\frac{\partial}{\partial\theta} \ell(\theta; X_i)\right] = 0$$

*and*

$$\text{E}_{\theta}\left\{\left[\frac{\partial}{\partial\theta}\log f(X_i; \theta)\right]^2\right\} = -\text{E}_{\theta}\left[\frac{\partial^2}{\partial\theta^2}\ell(\theta; X_i)\right]$$

*for $\theta\in \mathcal O$.*


The first derivative of the log-likelihood function

$$\frac{\partial}{\partial\theta} \ell(\theta; X_i)$$

is called the **score function**. The expectation of the squared score

$$\mathcal I(\theta) = \text{E}_{\theta}\left\{\left[\frac{\partial}{\partial\theta}\log f(X_i; \theta)\right]^2\right\}$$

is called the **Fisher information** and the second equation in the Lemma 1 the **Fisher Information Equality**. Now, we prove the lemma.

---
We start from the identity

$$1= \int f(x_i;\theta) dx_i $$

where integration over $\mathcal X$ is understood. By differentiating both sides of the equation with respect to $\theta$, and relying on **R6**, we obtain

$$0=\int \frac{\partial}{\partial \theta}f(x_i; \theta)dx_i  = \int \frac{\partial \ell(\theta; x_i)}{\partial \theta} f(x_i; \theta)dx_i = \text{E}_{\theta}\left[\frac{\partial \ell(\theta; X_i)}{\partial\theta}\right],$$

which proves the first result. (Notice that this result would also break down if **R1** is not met, unless the density $f(x_i;\theta)$ converges to zero at the end-points of its support for all $\theta\in\mathcal O$.)


As for the second result, we again rely on **R6** and we differentiate both sides once again, which gives

$$0 = \int \frac{\partial^2 \ell(\theta; x_i)}{\partial\theta^2} f(x_i;\theta)dx_i + \int\left(\frac{\partial \ell(\theta; x_i)}{\partial}\right)^2f(x_i;\theta)dx_i.$$

So,

$$\text{E}_{\theta}\left\{\left[\frac{\partial \ell(\theta; X_i)}{\partial\theta}\right]^2\right\} = -\text{E}_{\theta}\left[\frac{\partial^2 \ell(\theta; X_i)}{\partial\theta^2}\right].$$

<div style="text-align: right"> Q.E.D </div>
---

A straight forward corollary gives us the variance of the score function:

**Corollary 3** *Under the conditions of Lemma 1,*

$$\text{Var}_{\theta}\left[\frac{\partial \ell(\theta; X_i)}{\partial\theta}\right] = \mathcal I(\theta).$$



Lastly, we need three additional assumptions to obtain the asymptotic distribution of the MLE.

**Assumption (Regularity Conditions, cont)**
- **R7** *The density $f(x_i;\theta)$ is three times differentiable with respect to $\theta \in \mathcal O$.*
- **R8** *The Fisher information satisfies $0 < \mathcal I(\theta) < \infty$ on $\mathcal O$.*
- **R9** *For all $\theta$, there exists a constant $c$ and a function $M(x)$ (that may depend on $\theta$) such that*

    $$ \left \vert \frac{\partial^3 \log f(x_i; \theta)}{\partial \theta^3} \right \vert < M(x),$$  

    *with*  

    $$ \text{E}_{\theta_0}[ M(X_i) ] < \infty,$$  

    *for all $\theta_0 - c < \theta < \theta_0+c$ and almost all $x_i$.*


As before, **R7** implies its weaker counterpart **R5**. **R7** provides the conditions to Taylor-expand the first derivative of the log-likelihood around $\theta\_0$, truncate it after the first term, and ensure that the mean-value form (sometimes called the Lagrange form) of the remainder term is well-defined (See Theorem A1 and Corollary A1). As it will be clear from the statement of the following theorem, **R8** is assumed as its violation would lead to limiting distributions that have infinite variance or that degenerate to a single point (and if the second derivative is unbounded, the third derivative would not exist). The importance of condition **R9** will become apparent in the prove of our final theorem.

**Theorem 4** *Under conditions **R0** to **R9**, any consistent sequence $\hat\theta_n = \hat\theta_n(\mathbf X_n)$ of roots of the likelihood equation satisfies*

$$\sqrt n (\hat\theta_n - \theta) \rightarrow W$$

*in distribution as $n\rightarrow\infty$,where $W \sim \text{Normal}(0, \mathcal I(\theta)^{-1})$, i.e.,  $\text{E}[W]= 0$ and $\text{Var}[W] = \mathcal I(\theta)^{-1}$.*

---
By **R7**, we can rely on Taylor's Reminder Theorem to expand the score function around $\theta_0$ for any fixed $\mathbf x_n \in \mathcal X$ as

$$\ell'(\hat\theta_n;\mathbf x_n) = \ell'(\theta_0;\mathbf x_n) + (\hat\theta_n - \theta_0)\ell''(\theta_0;\mathbf x_n) + \frac{1}{2}(\hat\theta_n - \theta_0)^2\ell'''(\tilde\theta_n;\mathbf x_n)$$

where $\tilde\theta_n$ lies between $\theta_0$. But $\hat\theta_n$ is, by assumption, a root of the likelihood equation. Thus, the left-hand side is zero, and by rearranging terms, we obtain

$$\hat\theta_n - \theta_0 = \frac{\ell'(\theta_0;\mathbf x_n)}{-\ell''(\theta_0;\mathbf x_n) - \frac{1}{2}(\hat\theta_n - \theta_0)\ell'''(\tilde\theta_n;\mathbf x_n)}$$

and, multiplying both sides by $\sqrt{n}$,

$$\sqrt n (\hat\theta_n - \theta_0) = \frac{\frac{1}{\sqrt n}\ell'(\theta_0;\mathbf x_n)}{-\frac{1}{n} \ell''(\theta_0;\mathbf x_n) - \frac{1}{2n}(\hat\theta_n - \theta_0)\ell'''(\tilde\theta_n;\mathbf x_n)}.$$

We will show that

1. $n^{-1/2}\ell'(\theta_0;\mathbf X_n) \rightarrow V$ in distribution, where $V\sim \text{Normal}(0, \mathcal I(\theta_0))$;
2. $-n^{-1} \ell''(\theta_0;\mathbf X_n) \rightarrow \mathcal I(\theta_0)$ in probability; and
3. $\frac{1}{2n}(\hat\theta_n - \theta_0)\ell'''(\tilde \theta_n;\mathbf X_n)\rightarrow 0$ in probability.

We will prove **1.** to **3.** in order. Let $Y_i = \ell'(\theta_0; X_i) = (\partial/\partial \theta)\log f(X_i;\theta)\vert_{\theta = \theta_0}$, i.e., the score evaluated at $X_i$ and $\theta = \theta_0$. Then, from Lemma 1 and Corollary 3, $\text{E}\_{\theta_0}[Y_i] = 0$ and $\text{Var}\_{\theta_0}[Y_i] = \mathcal I(\theta_0),$ where $0<\mathcal I(\theta_0)<\infty$ by **R8**. Hence,

$$
\begin{aligned}
\frac{1}{\sqrt n}\ell'(\theta_0;\mathbf X_n) &= \frac{1}{\sqrt n}\sum_{i=1}^n Y_i=\sqrt{\mathcal I(\theta_0)} \left(\frac{n^{-1} \sum_{i=1}^nY_i}{\sqrt{\mathcal I(\theta_0)}/\sqrt{n}}\right) \\
&\longrightarrow \sqrt{\mathcal I(\theta_0)}Z
\end{aligned}
$$

in distribution as $n\rightarrow \infty$ by the Central Limit Theorem, where $Z \sim \text{Normal}(0,1)$. Hence, $n^{-1/2}\ell'(\theta_0;\mathbf X_n) \rightarrow V$ in distribution, where $V\sim \text{Normal}(0, \mathcal I(\theta_0))$, which proves **1.**

Next, it follows from the Weak Law of Large Numbers that

$$-\frac{1}{n}\ell''(\theta_0; \mathbf X_n) = -\frac{1}{n}\sum_{i=1}^n \ell''(\theta_0; X_i)\rightarrow -\text{E}_{\theta_0}[\ell''(\theta_0;X_i)]$$

in probability. Yet, by the Fisher Information Equality (see Lemma 1),

$$ -\text{E}_{\theta_0}[\ell''(\theta_0;X_i)] = \mathcal I(\theta_0)$$

which proves **2.**

Lastly, to show that $\frac{1}{2n}(\hat\theta_n - \theta_0)\ell'''(\tilde \theta_n;\mathbf X_n)\rightarrow 0$ in probability, it will suffice to show that $\frac{1}{n}\ell'''(\tilde\theta_n, \mathbf X_n)$ is bounded with probability tending to one as $n\rightarrow\infty$ since $\hat\theta_n$ is, by assumption, a consistent sequence such that $(\hat\theta_n - \theta_0) \rightarrow 0$ in probability. By **R9**, there exists $c>0$ and a function $M: \mathbb R \rightarrow\mathbb R$, such that for $\vert\tilde\theta_n - \theta_0\vert < c$

$$\left\vert\frac{1}{n}\ell'''(\tilde\theta_n, \mathbf X_n)\right\vert= \left\vert\frac{1}{n}\sum_{i=1}^n \ell'''(\tilde\theta_n; X_i)\right\vert < \frac{1}{n}\sum_{i=1}^n M(X_i)$$

almost surely and where

$$\text{E}_{\theta_0}[M(X_i)] < \infty.$$

But $(\hat\theta_n - \theta_0) \rightarrow 0$ in probability implies that $(\tilde\theta_n - \theta_0) \rightarrow 0$ in probability as well, since  $\tilde\theta_n$ lies between $\hat\theta_n$ and $\theta_0$ for all $n$. Hence, for an arbitrary $\epsilon > 0$, we can find $N_1$ such that $n> N_1$ implies that

$$P_{\theta_0}[\vert\tilde\theta_n - \theta_0\vert > c] < \epsilon,$$

which shows that for sufficiently large $n$ $\text{E}\_{\theta_0}[M(X_i)]$ is bounded. Further, by the boundedness of the expectation, the Weak Law of Large Numbers applies and $n^{-1}\sum_{i=1}^n M(X_i) \rightarrow \text{E}\_{\theta_0}[M(X_i)]$ with probability tending to one. Hence, for some $\gamma > 0$, we can find a $N_2$ such that $n> N_2$ implies

$$P_{\theta_0}\left[\left\vert\frac{1}{n}\sum_{i=1}^n M(X_i)\right\vert  > \text{E}_{\theta_0}[M(X_i)]+ \gamma\right] < \epsilon.$$

It follows that for $n> \max\\{N_1, N_2\\}$,

$$P_{\theta_0}\left[ \left\vert\frac{1}{n}\ell'''(\tilde\theta_n; \mathbf X_n)\right\vert >  \text{E}_{\theta_0}[M(X_i)]+ \gamma\right] < \epsilon.$$

As $\epsilon$ was arbitrary, this shows that $n^{-1}\ell'''(\tilde\theta_n)$ will be bounded with probability converging to one as $n\rightarrow\infty$, which proves **3.**

**1.** to **3.** show that the numerator of

$$\sqrt n (\hat\theta_n - \theta_0) = \frac{\frac{1}{\sqrt n}\ell'(\theta_0;\mathbf x_n)}{-\frac{1}{n} \ell''(\theta_0;\mathbf x_n) - \frac{1}{2n}(\hat\theta_n - \theta_0)\ell'''(\tilde\theta_n;\mathbf x_n)}$$

converges in distribution to $V \sim \text{Normal}(0, \mathcal I(\theta_0))$ while the denominator converges in probability to $\mathcal I(\theta_0)^{-1}$. By Slutsky's Theorem, it follows that

$$\sqrt n (\hat\theta_n - \theta_0) \rightarrow W $$

in distribution, where $W$ is Normally distributed with mean

$$\text{E}\_{\theta_0}[W] =\mathcal I(\theta_0)^{-1} \text{E}\_{\theta_0}[V] = 0$$

and variance

$$\text{Var}\_{\theta_0}[W]= \mathcal I(\theta_0)^{-2}\text{Var}\_{\theta_0}[V] = \mathcal I(\theta_0)^{-1}.$$

We are done.

<div style="text-align: right"> Q.E.D </div>
---

In sum, when we assume conditions **A1** to **A2** and **R0** to **R9**, the MLE will be consistent and has a Normal limiting distribution. In fact, the variance of the MLE will coincide with the **Rao-Cramer lower bound**, which is the smallest attainable variance among all asymptotically unbiased estimators.

<br/>

<h3> Appendix </h3>

In deriving the asymptotic distribution of the MLE, we relied on Taylor's Theorem. Although the theorem is well-known and appears in introductory calculus books, it is restated in this appendix with proof.

**Theorem A1 (Taylor's Theorem with mean-value form remainder)** *Let $f: \mathbb R \rightarrow \mathbb R$ be a function which is $k+1$ times differentiable in the interior of the closed interval $I$ between $x_0$ and $x$, where $x,x_0 \in \mathbb R$. Let $f^{(k)}$ be continuous on $I$ and let $g:\mathbb R \rightarrow R$ be continuous on $I$ with non-vanishing derivative in the interior of $I$. Then,*

$$\begin{aligned}
f(x) &= f(x_0) + f'(x_0) (x-x_0) + \cdots + \frac{1}{k!}f^{(k)}(x_0)(x-x_0)^d + R_k(x,x_0) \\
&= \sum_{l=0}^k \frac{f^{(l)}(x_0)(x-x_0)^l}{l!} + R_k(x,x_0) \\
&= P_k(x) + R_k(x, x_0),
\end{aligned}$$

*where*

$$P_k(x) = \sum_{l=0}^k \frac{f^{(l)}(x_0)(x-x_0)^l}{l!}$$

*is called the **$k$th-order Taylor polynomial** of the function $f$ at point $x_0$ and $R_d(x,x_0)$ is called the **remainder term**, and where remainder satisfies*

$$R_k(x,x_0) = \frac{f^{(k+1})(\xi)}{k!}(x - \xi)^k \frac{g(x) - g(a)}{g'(\xi)},$$

*for some real number $\xi$ that lies b*etween $x$ and $x_0$.*


Before proving the theorem, we'll need Cauchy's Mean value Theorem (sometimes called the Generalized Mean Value Theorem).

**Lemma A1 (Cauchy's Mean Value Theorem; Generalized Mean Value Theorem)** *If $f$ and $g$ are continuous real-valued functions on $I =[a, b]$ which are differentiable in the interior of $I$, then there exists a point $\xi \in (a, b)$ such that*

$$[f(b) - f(a)]g'(\xi) = [g(b) - g(a)*) f'(\xi).$$

---
For $x\in I$, define

$$h(x) = [f(b) - f(a)]g(x) - [g(b) - g(a)]f(x).$$

As the linear combination of differentiable functions is also differentiable, $h$ is differentiable in the interior of $I$. 	To prove the theorem, we have to show that $h'(\xi) = 0$ for some $\xi\in (a, b)$.

We first notice that

$$h(a) = f(b)g(a) - f(a)g(b) = h(b).$$

Hence, if $h$ is constant, $h' = 0$ on $I$ and the theorem holds for any any value $x\in (a, b)$. So, assume otherwise.


By the continuity of $h$, $h(a)$ has to be either a minimum or maximum of $h$ on $I$. First, suppose that $h(a) = h(b)$ is a minimum. Then, by the Extreme Value Theorem, there exists a point $\xi \in (a, b)$ at which $h$ attains a maximum and $h'(\xi) = 0$. Similarly, if $h(a)$ is a maximum, choose $\xi$ to be the point at which $h$ attains its minimum. This completes the proof.
<div style="text-align: right"> Q.E.D </div>
---

Notice what is called "the'' Mean Value Theorem is a special case of the Cauchy version with $g(x) = x$. Next, we prove the Theorem.

---
Without loss of generality assume $x > x_0$. For $u \in (x_0, x)$, define

$$Q_k(u) = f(u) + f'(u)(x - u) + \frac{f''(u)(x- u)^2}{2}+ \cdots + \frac{f^{(k)}(u)(x-u)^k}{k!}.$$

Notice that, by assumption on the existence of the $k+1$th derivative of $f$, $Q_k(u)$ is differentiable in the interior of $I$. Hence, by the Generalized Mean Value Theorem, there exists $\xi \in (x_0, x)$ such that

$$\frac{Q_k'(\xi)}{g'(\xi)} = \frac{Q_k(x) - Q_k(x_0)}{g(x) - g(x_0)}.$$

Since $Q_k(x_0) = P_k(x)$ and $Q_k(x) = f(x)$, $Q_k(x) - Q_k(x_0) = R_k(x, x_0)$. But

$$\begin{aligned}
\frac{d}{du}Q_k(u) &= f'(u) + \Bigg[f''(u)(x - u) - f'(u)\Bigg] + \left[\frac{f'''(u)(x - u)^2}{2} - f''(u)(x - u)\right] +\\
&\quad \left[\frac{f^{(4)}(u)(x - u)^4}{4!} - \frac{f^{(3)}(u)(x - u)^{2}}{2!}\right] +\\
&\quad  \cdots + \\
& \quad \left[\frac{f^{(k+1)}(u)(x - u)^k}{k!} - \frac{f^{(k)}(u)(x - u)^{k-1}}{(k - 1)!}\right] \\
&= \frac{f^{(k+1)}(u)(x - u)^k}{k!}.
\end{aligned}$$

Therefore,

$$[g'(\xi)]^{-1} \frac{f^{(k+1)}(\xi)(x - \xi)^k}{k!} = \frac{R_{k}(x, x_0)}{g(x) - g(x_0)}$$

and rearranging yields

$$R_k(x,x_0) = \frac{f^{(k+1)}(\xi)(x - \xi)^k}{k!}\frac{g(x) - g(x_0)}{g'(\xi)}$$

as desired.
<div style="text-align: right"> Q.E.D </div>
---

**Corollary A1 (Lagrange Form of Remainder)** *Under the assumptions of Theorem A1,*

$$R_k(x, x_0) = \frac{f^{(k+1)}(\xi)(x - x_0)^{k+1}}{(k+ 1)!}.$$

---
Take $g(u) = (x - u)^{k+1}$. Then

$$\frac{g(x_0) - g(x)}{g'(\xi)} = \frac{(x- x_0)^{k + 1}}{(k + 1)(x - \xi)^k}$$

and the result follows.
<div style="text-align: right"> Q.E.D </div>
---

Another theorem that was used in the derivation of the asymptotic distribution of the MLE was Slutsky's Theorem. This theorem appears in various forms in probability textbooks. Here, we'll state and prove the version that was used in the last proof.


We start with a lemma that will not be proven.

**Lemma A2** *Let $\\{X\_n\\}\_{n\ge 1}$ be a sequence of real-valued random variables. $X\_n \rightarrow X$ in distribution if and only if $\lim\_{n\rightarrow\infty}\text{E}[g(x\_n)] = \text{E}[g(X)]$ for all bounded Lipschitz continuous functions $g$.*

In some text books, the following theorem is introduced as "the'' Slutsky's Theorem.

**Theorem (Slutsky's Theorem)** *Let $\\{X\_n\\}\_{n\ge1}$ and $\\{Y\_n\\}\_{n\ge 1}$ be two sequences of$\mathbb R^d$ valued random variables, with $X\_n \rightarrow X$ in distribution and $\|X\_n - Y\_n\| \rightarrow 0$ in probability. Then $Y\_n \rightarrow X$ in distribution.*

---
By Lemma A2, it suffice to show that $\lim_{n\rightarrow\infty}\text{E}[f(Y_n)] = \text{E}[f(X)]$ for all Lipschitz continuous, bounded, $f$. Let $f$ have these properties. Then, there exist positive constants $k$ and $\alpha$ such that $\vert f(x) - f(y)\vert\le k \|x - y\|$ and $\sup_{x}\vert f(x)\vert < \alpha$. Hence

$$
\begin{aligned}
\lim_{n\rightarrow\infty}\vert\text{E}[f(X_n) - f(Y_n)]\vert & \le \lim_{n\rightarrow\infty}\text{E}\Big[\Big\vert f(X_n) - f(Y_n)\Big\vert\Big] & \text{(Jensen's Inequality)} \\
&=\lim_{n\rightarrow\infty}\text{E}\Big[\Big\vert f(X_n) - f(Y_n)1_{\{\|X_n - Y_n\| \le \epsilon\}} & \\
&\qquad \qquad \qquad + f(X_n) - f(Y_n)1_{\{\|X_n - Y_n\| > \epsilon\}}\Big\vert\Big]&\\
&\le\lim_{n\rightarrow\infty}\text{E}\Big[\Big\vert f(X_n) - f(Y_n)1_{\{\|X_n - Y_n\| \le \epsilon\}}\Big\vert & \text{(Triangle Inequality)}\\
&\qquad + \Big\vert f(X_n) - f(Y_n)1_{\{\|X_n - Y_n\| > \epsilon\}}\Big\vert\Big]&\\
&\le k\lim_{n\rightarrow\infty}\text{E}\Big[\|X_n - Y_n\|1_{\{\|X_n - Y_n\| \le \epsilon\}}\Big\vert & \text{(Lipschitz cont of $f.$)}\\
&\qquad  + \lim_{n\rightarrow\infty}\text{E}\Big[\Big\vert f(X_n) - f(Y_n)1_{\{\|X_n - Y_n\| > \epsilon\}}\Big\vert\Big]&\\
&\le k\epsilon + \lim_{n\rightarrow\infty}\text{E}\Big[\Big\vert f(X_n) - f(Y_n)1_{\{\|X_n - Y_n\| > \epsilon\}}\Big\vert\Big]&\\
&\le k\epsilon + \lim_{n\rightarrow\infty}2\alpha\text{E}\Big[1_{\{\|X_n - Y_n\| > \epsilon\}}\Big]&\text{(Boundedness of $f$)}\\
&= k\epsilon + 2\alpha\lim_{n\rightarrow\infty}P[\|X_n - Y_n\| > \epsilon]&\\
&\le\epsilon( k + 2\alpha),
\end{aligned}
$$

where the last step follows from $X_n \rightarrow Y_n$ in probability. As $\epsilon$ is arbitrary, it follows that

$$\lim_{n\rightarrow\infty}\text{E}[f(Y_n)] = \lim_{n\rightarrow\infty}\text{E}[f(X_n)] = \text{E}[f(X)]$$

for all Lipschitz continuous, bounded functions $f$. We are done.
<div style="text-align: right"> Q.E.D </div>
---

Other textbooks introduce Slutsky's Theorem in a different form, which is more useful for many applications. As we have already proven the first version, proving the other version will be straightforward.

**Theorem A3 (Slutsky's Theorem)** *Let $X_n$ converge in distribution to $X$ and let $Y_n$ converge in probability to a constant $c$.* Then,

$$\begin{aligned}
X_n + Y_n &\rightarrow X + c \\
X_nY_n & \rightarrow cX \\
X_n/Y_n & \rightarrow X/c
\end{aligned}$$

*in distribution.*

---
We will first show that $(X_n , Y_n)$ converges in distribution to $(X, c)$. Let $f(x, y)$ be a bounded, Lipschitz continuous function and let $g: x\mapsto f(x, c)$. Clearly $g$ is bounded and Lipschitz as well. But as $X_n \rightarrow X$ in distribution,
$$\lim_{n\rightarrow\infty}\text{E}[f(X_n, c)] = \lim_{n\rightarrow\infty}\text{E}[g(X_n)] = \text{E}[g(X)] = \text{E}[f(X, c)].$$

Hence $(X_n, c) \rightarrow (X, c)$ in distribution.

Next, consider $\|(X_n, Y_n) - (X_n, c)\|_1 = \vert X_n - X_n\vert + \vert Y_n - c\vert = \vert Y_n - c\vert$. As $Y_n \rightarrow c$ in probability, $\vert Y_n - c\vert \rightarrow 0$ in probability as well. So, we have established that
$\|(X_n, Y_n) - (X_n, c)\| \rightarrow 0$ in probability and $(X_n, c) \rightarrow (X, c)$ in distribution. By (first version of) Slutsky's Theorem from above, this implies that $(X_n, Y_n) \rightarrow (X, c)$ in distribution.

Lastly, we notice that the functions $h_1: (x, y) \mapsto x + y,$ $h_2:(x, y) \mapsto xy$, and $h_3: (x, y) \mapsto x/y$ are all continuous. Hence, by the Continuous Mapping Theorem, the desired results follow.
<div style="text-align: right"> Q.E.D </div>
---
