---
title : Seriously, What is a Random Variable?
date : 2019-09-07 02:44:02
author : baruuum
comments : true
---


So, I'm the TA for the Introduction to Statistics class for the incoming PhD cohort. So far, every statistics course that I've taken in sociology departments introduced a random variable as "a variable that is random." In the the class I am TAing for, on the other hand, the concept of a (real-valued) random variable was introduced quite rigorously, namely as a mapping from the sample space to the real line. I loved it, but the students seemed to be lost. After preparing [a document]({{ site.baseurl }}{% link _teaching/002_2019_MathRefresh.md %}) to help students to walk through the definitions and some examples, I thought that it was not quite complete. So, in this post, I cram together my rusty and limited knowledge of probability theory to summarize what a random variable is.

----

### Probability Spaces

The definition of a random variable starts necessarily from **probability spaces**. A probability space is the triple $(\Omega, \mathcal F, P)$ (Why is this called a space? Since mathematicians call many things with some kind of structure a space.) We might think of the probability space as a *model* of a random phenomenon. It is a model, because it is based on a set of assumptions and restrictions about the underlying structure we try to study. Let's discuss each element of this triple in turn.

1. $\Omega$ denotes the **sample space**. This is the set of all possible outcomes of an experiment. Notice that this entails the assumption that *we know* what the possible outcomes of the experiment are. For example, in the case of throwing a single coin, we have $\Omega = \\{H, T\\}$, where $H$ stands for Heads and $T$ for Tails. But it could also happen that the coin lands on its edge; and we are assuming that this outcome is not part of the random process we are studying (Of course, we could extend our sample space to $\Omega = \\{H, T, E\\}$, where $E$ stands for the coin landing on its edge.)

2. $\mathcal F$ is *a* (not *the*) collection of subsets of $\Omega$, which we call **events**. Notice that $\mathcal F$ is a set of sets (usually we use other terms, such as *collection* or *family*, to denote sets whose elements are sets as well). $\mathcal F$ is sometimes called the **event space**. We say that an event $A \in \mathcal F$ **occurs** if as a result of the experiment we observe some outcome $\omega\in A$. In other words, if we throw a coin and get the result $H\in \Omega$, and the set $\\{H\\}$ belongs to the collection $\mathcal F$, then the event $\\{H\\}$ has occurred (Notice that the trivial event, $\Omega = \\{H,T\\}$, has occurred as well as $H\in \\{H, T\\}$.)

    But here it gets a little bit complicated. While it is tempting to define arbitrary subsets of $\Omega$ as events, it turns out that this can lead to contradictions (such as the [Banach-Tarski Paradox](https://en.wikipedia.org/wiki/Banach%E2%80%93Tarski_paradox)). So, from the outset, we restrict the collection of events to subsets which behave in a reasonable way. The restriction we impose on the event space, $\mathcal F$, is that it forms a **$\sigma$-algebra** on $\Omega$. While the term sounds scary, it is just a structure that tells us what subsets of $\Omega$ should be counted as events. The conditions are the following:
    
    1. $\Omega \in \mathcal F$. 
        + The sample space is an event, sometimes called the *sure event*.
    2. If $A \in \mathcal F$, then $A^C = \Omega\setminus A \in \mathcal F$. 
        + If a $A$ is an event, then its relative complement is also an event. This basically means that if an event can happen, there must be another event that represents the scenario in which it does not happen.
        + As $\Omega \in \mathcal F$ by Condition 1., this implies that $\varnothing \in \mathcal F$ as well.
    3. Let $(A\_n)\_{n=1}^\infty$ be a sequence of events, i.e., $A\_n \in \mathcal F$ for all $n\ge 1$. Then $\bigcup\_{n=1}^\infty A\_n \in \mathcal F$. 
        + The countable union of events is also an event. Heuristically, this means that for every (countable) sequence of events, "at least one of the events in the sequence happens" is also an event.
        + Here, *countable* means that we can *enumerate* the elements in our sequence as $\\{A\_1,A\_2, A\_3,...\\}$. For example, consider the set of natural numbers, $\mathbb N$. If you had infinite time (where time proceeds in discrete steps), you could enumerate all natural numbers, starting with $1$, then $2$, and so on. So the set of natural numbers, $\mathbb N$, is *enumerable* but infinitely large. We say in these situations that $\mathbb N$ is *countably infinite*. On the other hand, there are sets that are "too large" for such enumerations. The set of real numbers in the interval $[1,a)$, $a> 1$ would be such an example. Regardless how you choose your next step after the number $1$, say $1 + \epsilon$, where $0<\epsilon < a - 1,$ there will be, again, infinitely many real numbers between $1$ and $1+\epsilon$. In fact, there will be more real numbers in the interval $[1, 1+\epsilon)$, then there are natural numbers in $\mathbb N$; it is a different kind of infinity, which we call *uncountable*.
    
    There are not many meaningful countable unions for our coin-tossing experiment. But just for the sake of showing an example, let $\\{H\\}\in \mathcal F$. Then by Condition 2. $\\{T\\} = \\{H\\}^C \in \mathcal F$ as well. Next, consider the sequence $(A\_n)\_{n=1}^\infty$, where $A\_n = \\{H\\}$ if $n$ is even and $A\_n = \\{T\\}$ otherwise. Condition 3. requires that $A = \bigcup\_{n=1}^\infty A\_n$ is a set in $\mathcal F$. For the current example, this is vacuously true as $A = \\{H, T\\} = \Omega$, which must be a member of $\mathcal F$ by condition 1. Hence, $\sigma$-algebra of our coin tossing experiment would look like $\mathcal F = \\{\varnothing, \\{H\\}, \\{T\\}, \\{H, T\\}\\}$. Notice that we could also started by saying that $\\{H\\}$ is *not* in $\mathcal F$. Then, we would get the trivial $\sigma$-algebra $\mathcal F = \\{\varnothing, \Omega\\}$.
    
    As for another example, suppose we draw a random natural number. In this case, the sample space is $\Omega = \mathbb N$. Let's treat all subsets of $\Omega$ of the form $\\{i\\}, i \in \mathbb N$, as events. Condition 3. would then require that the union $O = \bigcup\_{n=1}^\infty O\_n$, where $O\_n = \\{\text{$n$th smallest odd natural number}\\}$, must be a set in $\mathcal F$, since it is a countable union of the elementary events $\\{i\\}, i\in \mathbb N$. Further, by condition 2., we have that the set of even numbers must be an event as well. Hence, if the elementary events of drawing any natural number are treated as events, drawing an odd number or an even number have to be events as well by the restrictions we impose on the event space. 
    
3. The pair $(\Omega, \mathcal F)$ forms a **measurable space**. Again, the words *measurable* and *space* sound scary, but it is just a set of objects together with a $\sigma$-algebra on that set, where the important point here is that $\mathcal F$ is a $\sigma$-algebra. The sets in $\mathcal F$ are called *measurable sets* (hence, events are measurable subsets of the sample space); these are the subsets of $\Omega$ which we want and can study/measure. 

    A **measure** is a set function from the $\sigma$-algebra of measurable sets to a non-negative real number or positive infinity, i.e.,

    $$\mu : \mathcal F \longrightarrow [0,\infty].$$
    
    Intuitively, a measure tells us "how large," so to speak, a subset of $\Omega$ is. What "how large" means depends on the measure we use. For example, it could mean "the number of elements" in a subset or it could mean the "length" of an interval of real numbers. When a measurable space is equipped with a measure, we call the triple $(\Omega, \mathcal F, \mu)$ a **measure space**. 

    To tell us the "size" of a set, it seems reasonable to ask a measure to satisfy some basic assumptions. For example, it is natural to ask for a function that measures the size of a set to be zero for the empty set; we might also require that it should not be negative; and for a *disjoint* collection of sets, by which we mean that the intersection of any pair of sets in the collection is empty, the measure of their union should be equal to their sum. For example, if we want to measure the size of set that consists of two intervals that don't overlap, it sounds reasonable to require that the size of the set should be equal to the sum of the lengths of the intervals. More formally, we can express these conditions as follows:
    
    1. $\mu(\varnothing) = 0$ (the measure of the empty set is zero).
    2. For all $A\in \mathcal F, \mu(A) \ge 0$ (the measure of a measurable set is non-negative).
    3. Let $(A\_n)\_{n = 1}^\infty$ be a sequence of *pairwise disjoint* measurable sets, i.e., $A\_n \cap A\_{n'} = \varnothing$ for all $n \ne n'$. Then $\mu(\bigcup\_{n=1}^\infty A\_n) = \sum\_{n=1}^\infty \mu(A\_n)$ (this is often called countable additivity).
        
4. Now, the last element of our probability space, $P$, is a **probability measure**. As it is a measure, it must satisfy the three conditions above. Further, for a probability measure, we require in addition that

    $$P(\Omega) = 1.$$

    In other words, the probability measure tells us the "size" of a subset of $\Omega$, with the requirement that the "size" of the set of all possible outcomes of an experiment is equal to one. This clarifies an important point that is often forgotten by students: probabilities are *not* functions from $\Omega$ to $[0,1]$. Probabilities are assigned to events, which are (measurable) subsets of $\Omega$. When flipping a coin, the probability of observing a Heads outcome would be $P(\\{H\\})$ not $P(H)$, although by writing $P(H)$ it is usually understood that we are referring to $\\{H\\}$. 


This, I hope, explains what a probability space is. To discuss random variables, we need one more piece: the definition of a **measurable function**.

Consider two measurable spaces, $(G, \mathcal G)$ and $(H, \mathcal H)$. A function 

$$f: G \longrightarrow H$$ 

is called **measurable** if for every $A \in \mathcal H$, the **preimage** under $f$, written as $f^{-1}(A)$, is $\mathcal G$-measurable. By preimage of $A$ under $f$ we mean the set of elements in $G$ that are mapped to $A$ by the function $f$, i.e., $f^{-1}(A) = \\{x \in G: f(x) \in A\\}$.  Notice that $A$ is a subset of $H$, while $f^{-1}(A)$ is a subset of $G$. By $\mathcal G$-measurable, we mean that $f^{-1}(A) \in \mathcal G$. Sometimes we write $f: (G, \mathcal G) \longrightarrow (H, \mathcal H)$ to emphasize the $\sigma$-algebras of measurable sets in both the domain and co-domain.

This definition says nothing about whether the image of a set $A \subseteq G$ under $f$ is measurable in $(H, \mathcal H)$. So, even if $f(A) \notin \mathcal H$, this will not affect the measurability of $f$. This sounds weird at first. But we'll see why measurable functions are defined in this way (at least this is one good reason that I can think of.) Namely, if our domain is equipped with a measure and the preimage of all measurable sets in the target space is measurable, then there is a natural way to define a new measure on the target space through the function $f$. 

This makes us ready to discuss what a random variable is.

### Random Variables and Distributions

Finally(!), a **random variable** $X$ is a measurable function that maps $(\Omega, \mathcal F)$ to another measurable space $(E, \mathcal E)$. In short,

$$ X: (\Omega, \mathcal F) \longrightarrow (E, \mathcal E).$$

Sometimes it is said that $X$ is a $(E, \mathcal E)$-valued random variable to emphasize the structure of the target space. By definition of measurable functions, the preimage of $A\in \mathcal E$ under $X$ is an event, i.e., $X^{-1}(A) \in \mathcal F$. This is important because we have already a well-defined probability measure on $(\Omega, \mathcal F)$, namely $P$. 

As we can assign probabilities to all sets in $\mathcal F$, this opens up the door to define probabilities to (measurable) subsets of the co-domain, $E$. We do this by equipping the target space $(E, \mathcal E)$ with a new probability measure $P\_X$ that is defined as

$$P_X(A) = P(X^{-1}(A)) = P(\{\omega \in \Omega: X(\omega) \in A\})$$

for all $A\in \mathcal E$ (notice the subscript on $P\_X$). In other words, for a set $A\in \mathcal E$, which is a set in the co-domain of the random variable $X$, we define the probability of $A$ as the probability of the set constructed of all $\omega \in \Omega$ such that $X(\omega) \in A$. As $X$ is by assumption a measurable function, $f^{-1}(A)$ is $\mathcal F$-measurable for all $A\in\mathcal E$, and thus $P_X$ is well-defined. 

A measure so constructed is often called a **pushforward** of the probability measure $P$ by $X$ and is sometimes denoted by $X\_\*\mu$. In the case where the domain of the function $X$ is a probability space, it is called a **distribution measure** or simply the **distribution** or **law** of the random variable $X$. 

Once we have defined a probability measure directly on $(E, \mathcal E)$, we can "forget" about the underlying probability space $(\Omega, \mathcal F, P)$ and work with the space "induced" by the random variable: $(E, \mathcal E, P\_X)$. The only thing we have to check is that $(E, \mathcal E, P\_X)$ is indeed a probability space. Since $(E, \mathcal E)$ is by definition a measurable space, it suffices to show that $P\_X$ is a probability measure on $(E, \mathcal E)$. This follows from the measurability of $X$ and the fact that preimages preserve set operations. Namely,

1. $P_X(\varnothing) = P(X^{-1}(\varnothing)) = P(\varnothing) = 0$ (measure of empty set is zero)
2. $P_X(E) = P(X^{-1}(E)) = P(\Omega) = 1$ (measure of universal set is one).
1. $P_X(A) = P(X^{-1}(A)) \ge 0$ for all $A\in \mathcal E$, as $P$ is non-negative on $\mathcal F$ (non-negativity).
3. As for countable additivity, let $(A_n)\_{n= 1}^\infty$ be a disjoint sequence in $\mathcal E$. Then,
$$\begin{aligned}
 P_X\left(\bigcup_{n= 1}^\infty A_n\right) &= P\left(X^{-1}\left[\bigcup_{n= 1}^\infty A_n\right]\right)= P\left(\bigcup_{n= 1}^\infty X^{-1}(A_n)\right) = \sum_{n= 1}^\infty P\Big(X^{-1}(A_n)\Big) \\
 & = \sum_{n= 1}^\infty P_X(A_n).
 \end{aligned}
$$ 
    
Thus, $P_X$ is a valid probability measure on $(E, \mathcal E)$ as desired.

Next, a **real-valued random variable** is a random variable where the co-domain is the real line.

$$ X: \Omega \longrightarrow \mathbb R.$$

Of course, we have to specify the $\sigma$-algebra of measurable sets in our co-domain. It is common to use the **Borel $\sigma$-algebra** on $\mathbb R$ for this purpose, denoted by $\mathcal B(\mathbb R)$ or simply $\mathcal B$. The Borel $\sigma$-algebera is the $\sigma$-algebra generated by the open intervals in $\mathbb R$. *Generated* here means that the Borel $\sigma$-algebra is the smallest $\sigma$-algebra that contains the open intervals in $\mathbb R$, where *smallest* means that it is contained in all $\sigma$-algebras that also contain the open intervals. For example, as the power set of $\mathbb R$ is a $\sigma$-algebra on $\mathbb R$, and since it contains *all* subsets of $\mathbb R$, it contains the open intervals as well. Hence, $\mathcal B(\mathbb R) \subseteq \mathcal P(\mathbb R)$ (i.e., the Borel $\sigma$-algebra is contained in the power set) and the former is *smaller* than the latter.

We have introduced the Borel $\sigma$-algebra as being generated by the open intervals in $\mathbb R$. But, as complements and countable unions have to be elements of $\mathcal B$ as well (it is a $\sigma$-algebra after all), it contains also all closed intervals and half-open intervals. This is because closed or half-open intervals can be expressed as a countable intersection of open intervals. Notice that this implies that all intervals of the form $(-\infty, x], x\in \mathbb R,$ will be sets in $\mathcal B$. Further, $X$ will be a random variable only if the preimage of $(-\infty, x], x\in \mathbb R,$ under $X$ are events, since $X$ must be measurable.

Probabilities assigned to intervals of the form $(-\infty, x]$ correspond to the **cumulative distribution function (CDF)** of the random variable $X$, 

$$F_X : \mathbb R \longrightarrow [0,1],$$

where 

$$F_X(x) = P_X((-\infty, x]) = P(\{\omega\in \Omega: X(\omega) \le x\}), \qquad x\in\mathbb R.$$

Further, it turns out that $F\_X$ uniquely defines $P\_X$. In other words, if we know $F\_X$, we know $P\_X$ and vice versa. The proof of this is quite technical but can be found in any standard probability textbook. The important point is that $F\_X$ exists whenever the distribution measure $P\_X$ exists and that it gives us a means to, again, "forget" about $P\_X$ and work with a function of real numbers, $F_X(\cdot)$.

----

Going from here to probability density functions needs another whole set of concepts and machinery (such as absolute continuity of measures and Radon-Nikodym derivatives). So, I rather stop here. But, before concluding, one might ask why we need random variables at all, given random experiments can be described perfectly with a probability space? 

I have no good answer except that it is just super convenient. Expressing random phenomena with $(\Omega, \mathcal F, P)$ is really tedious once these phenomena get complex. For example, I have no idea how one would prove any of the fundamental theorems in statistics, such as the Law of Large Numbers or the Central Limit Theorem, without using random variables. So, while all these definitions are important to work with random phenomena in a mathematically rigorous way, the ability to forget about the probability space and deal with random variables make life just much, much easier.
