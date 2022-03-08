---
title : Sampling Integers with RcppArmadillo
date : 2019-06-24 02:47:36
author : baruuum
comments : true
---

**DISCLAIMER**: _This is a blog post from my graduate school years. The blog is no longer maintained and the material might include typos._

To my surprise, it turns out that sampling integers with unequal probabilities using the `RcppArmadillo` the package is more subtle than it appeared at first. To explain the problem, consider the following situation. We have a length-$K$ vector, $\mathbf p = [p\_1,p\_2,...,p\_K]^\top$ which satisfies $\sum\_{j=1}^K p\_j = 1$. In a `C++` code that uses the [Armadillo library](http://arma.sourceforge.net/), we want to sample an integer from $\mathcal X = \\{1,2,...,K\\}$ according to these probabilities, i.e., sample the integer $1$ with probability $p\_1$, $2$ with probability $p\_2$ and so on. In addition, after sampling, we need the vector $\mathbf p$ for some further calculations and the end result of the program is to be exported into `R`.

While the sampling function of the `RcppArmadillo` package is super convenient, it has a small problem in these situations: namely, the function modifies the probability vector, $\mathbf p$,in its calculations. 

The following code might illustrate the point. First we load the `RcppArmadillo` package:


{% highlight r %}
library(RcppArmadillo)
{% endhighlight %}

Our code for sampling looks like the following. First, we include the header for the sampling function and tell `Rcpp` that we the code depends on the `Armadillo` library. It suffices to include the `RcppArmadilloExtensions/sample.h` header as it will include `RcppArmadillo.h` automatically.


{% highlight cpp %}
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
{% endhighlight %}

Next, we write our sampling function:


{% highlight cpp %}
//[[Rcpp::export]]
arma::uword sample_arma(arma::vec & pvec) {

    // count number of elements in pvec
    arma::uword K = pvec.n_elem;
    // create integers 0 to (K-1) to sample from
    arma::uvec opts = arma::linspace<arma::uvec>(0L, K - 1L, K);
    // sample integer
    return arma::conv_to<arma::uword>::from(
            Rcpp::RcppArmadillo::sample(opts, 1L, true, pvec)
    );

}
{% endhighlight %}
Notice that we are passing the `pvec` object by reference to the function in order to avoid making a deep copy. Second, as the `Rcpp::RcppArmadillo::sample` function returns an `unsigned int` object, we have to explicitly convert it into an `uword` object (which corresponds to an unsigned integer).

Now, when testing the code we observe the following:


{% highlight r %}
set.seed(123)
# no of integers to sample from
k = 4L
# prob vector
p = runif(k)
pvec = p / sum(p)
print(pvec)
{% endhighlight %}

```
## [1] 0.1214495 0.3329164 0.1727188 0.3729152
```

{% highlight r %}
# sample integers
res = sample_arma(pvec)

# look at results
print(res)
{% endhighlight %}

```
## [1] 0
```

{% highlight r %}
# look at prob vector
print(pvec)
{% endhighlight %}

```
## [1] 0.3729152 0.7058317 0.8785505 1.0000000
```
Observe that the `pvec` object has changed after the sampling is done: it now contains the cumulative probabilities. If we need to do some subsequent calculations with the `pvec` object after using it to sample integers, this might be problematic. We would face two options: add a calculate the finite differences to recover the original probability vector or pass the vector of probabilities by value. (Unfortunately, passing `pvec` by constant reference will not work). Both of these optinos might become expensive if the program needs to sample very often and `pvec` has many elements.

After having a look into the [source code](https://github.com/RcppCore/RcppArmadillo/blob/master/inst/include/RcppArmadilloExtensions/sample.h) of the function, I decided to write up a new function based on it. The logic of how the function works is quite straightforward. If we have only two possible outcomes, $1$ and $2$, and we want to sample $1$ with probability $p\_1$ and $2$ with probability $p\_2 = 1-p\_1$, we can just draw a random number, $u$, from the interval $(0,1)$ and return $1$ if $u < p\_1$ and $2$ otherwise. As the probability of $u$ being smaller or equal to $p\_1$ is exactly $p\_1$, we would be drawing $1$ according to the correct probability.

For the general case, we consider the length-$K$ vector of cumulative probabilities, $q = [q\_1, q\_2, ..., q\_K]$, where $q\_1 = p\_1$ and $q\_j = \sum\_{i \le j} p\_j$. Since  $q\_{j} - q\_{j-1} = \sum\_{i\le j}p\_j - \sum\_{i\le j-1} p\_j = p\_j$, we see that the length of the interval $(q\_j - q\_{j-1})$ is equal to the probability with which we want to sample the $j$th outcome. Hence, we can draw $U\sim \text{Uniform}(0,1)$ and choose outcome $j$ if $q\_{j-1} < U < q\_j$, where we define $q\_0$ to be equal to zero.

Now, there are two additional *tricks* that we might want to consider. First, after sampling $U$ from $(0,1)$, how would we find the interval into which $U$ falls? We might loop over all the intervals and stop when the desired interval is reached; yet, this can be inefficient if we want to sample multiple observations using the same probability vector, because some of the intervals, $(q\_j - q\_{j-1})$, are wider than the others. If we start the loop from the first interval and advance in order, it would be great if the first interval is the widest and, thus, has the largest probability of $U$ falling into it, the second interval the second largest, and so on. This could save us quite a lot of time, as for each draw, we likely would have to iterate over a smaller number of intervals. 

To incorporate this consideration, we might write our tentative sampling code as follows:


{% highlight cpp %}
//[[Rcpp::export]]
arma::uword samp_temp(const arma::vec &pvec) {
    
    // generate a vector of indices, so that 0 represents the largest 
    // element, 1 the second largest, and so on
    arma::uvec indx = arma::sort_index(pvec, "descend");
    
    // generate the q-vector, the vector of cumulative sums
    arma::vec qvec = arma::cumsum(arma::sort(pvec, "descend"));
    
    // draw randomly from (0,1)
    double u = arma::randu();
    
    // find interval into which u falls
    for (arma::uword k = 0; k < pvec.n_elem; ++k) {
        
        if (u <= qvec(k))
            return indx(k);
        
    }
}
{% endhighlight %}
Notice that once `u <= qvec(k)` is true, the program immediately terminates after returning the $k+1$th index (it is the $k+1$th index because indices in `C++` start from zero). Thus, if the loop does not terminate in the first iteration (`k=0`), we know that `u > qvec(0)`; if it then terminates in the second iteration (`k=1`), this means that `qvec(0) < u < qvec(1)`, and so on. Also, the code will draw only a single integer; the extension to drawing `N` samples would be straightforward (although we'll not pursue this here).

This is all good, but these calculations can become quite unstable if the probabilities contained in `pvec` get very small. So, for better numerical stability, it is always good to operate in the log-space. We won't go into too much detail, but only note that there is a special function to calculate $\log(1 + x)$ in a numerically stable way. Thus, whenever you encounter a situation in which you have to add one to a real number and take the logarithm, never use `log(1 + x)` but instead use `log1p(x)`, which is provided in almost all programming languages (for example, try typing `?log1p` into your `R` console). 

So, why is this important? Suppose we are working with log-probabilities. We have $\log p\_i, \log p\_j \ll 0$, and we want to calculate $\log(p\_i + p\_j)$. As both $\log p\_i$ and $\log p\_j$ are very small, exponentiating them might lead to underflow. For example, if you have $\log p\_i = \log p\_j = -750$ and you exponentiate them in `R`, you'll get;

{% highlight r %}
lp1 = lp2 = -750
print(exp(lp1))
{% endhighlight %}

```
## [1] 0
```
Hence, `exp(lp1) + exp(lp2)` will evaluate to zero, and taking the logarithm evaluates to `-Inf`. 

{% highlight r %}
print(log(exp(lp1) + exp(lp2)))
{% endhighlight %}

```
## [1] -Inf
```
We know that this answer cannot be true: since both $p\_1$ and $p\_2$ are finite, their sum is finite as well, and $\log(p\_1 + p\_2) > -\infty$. 

A better way to calculate $\log(p\_i + p\_j)$ is the following. We rewrite the quantity as


$$\begin{aligned}
\log(p_i + p_j) & =\log(e^{\log p_i}+ e^{\log p_j}) = \log\left[ p_i\Big(1 + p_i^{-1}e^{\log p_j}\Big) \right] = \log p_i + \log\Big(1 + e^{\log p_j - \log p_i}\Big) \\
&=\log p_i + \text{log1p}\Big[\exp(\log p_j - \log p_i)\Big]
\end{aligned}$$


and use the `log1p` function to evaulate $\log(p\_1 + p\_2)$ as follows:

{% highlight r %}
print(lp1 + log1p(exp(lp2 - lp1)))
{% endhighlight %}

```
## [1] -749.3069
```
In `C++` we might code up a function to calculate the logarithm of the sum of two exponentials as follows:

{% highlight cpp %}
inline double log_add_exp(const double x, const double y) {

    double d = x - y;
    if (d > 0.0)
        return x + std::log1p(std::exp(-d));
    if (d <= 0.0)
        return y + std::log1p(std::exp(d));

}
{% endhighlight %}
We might, further, use this function to accumulate log-probabilities on the log-scale:

{% highlight cpp %}
template <typename T>
inline T log_accu_exp (const T& x) {

    // initialize with x
    T y(x);
  
    // accumulate probabilities on the log-scale
    typename T::iterator i = y.begin() + 1;
    for (; i < y.end(); i++) {
        *i = log_add_exp(*(i - 1), *i);
    }

    return y;
}
{% endhighlight %}
So, if we enter a vector $\mathbf x = [\log(x\_1),\log(x\_2),...,\log(x\_K)]^\top$ into `log_accu_exp`, the first element of the returned vector will be $\log(e^{\log(x\_1)}) = \log(x\_1)$, the second element $\log(e^{\log(x\_1)} + e^{\log(x\_2)}) = \log(x\_1 + x\_2)$, and so on.

Now, we are ready to write our second (tentative) sampling function. Here we use the fact that if $U\sim \text{Uniform}(0,1)$, then for a constant $c>0$, 



$$ \Pr[ -\log(U) \le c] = 1 - \Pr[U \le \exp(-c)] = 1 - e^{-c},$$



which shows that $Z = -\log(U)$ is distributed as a $\text{Exponential}(1)$ random variable. To sample the log of a $\text{Uniform(0,1)}$ distribution, we can therefore sample from the $\text{Exponential}(1)$ distribution and negate the result.


{% highlight cpp %}
//[[Rcpp::export]]
arma::uword lsamp_temp(const arma::vec &lpvec) {
    
    //check whether all elements are finite
    if (lpvec.has_inf())
        Rcpp::stop("log-probabilities have to be finite");

    // get indices of elements in descending order
    arma::uvec indx = arma::sort_index(lpvec, "descend");
    // logarithm of accumulated probabilities
    arma::vec lqvec = log_accu_exp(arma::vec(arma::sort(lpvec, "descend")));

    // sample log(unif(0,1))
    double r = -R::rexp(1.0);

    // sample integer
    for (arma::uword k = 0; k < lqvec.n_elem; ++k) {
        if (r <= lqvec(k))
            return indx(k);

    }

}
{% endhighlight %}

Lastly, we accomodate the function so that it can take log-probability vectors which are not properly normalized. Let $\mathbf p^{\ast} = [p\_1^{\ast}, ..., p\_k^{\ast}]^\top$ be a unnormalized probability vector such that $\mathbf p = s^{-1}\mathbf p^{\ast}$ and $s = \sum\_j p\_j^{\ast}$. To sample from $\mathbf p^{\ast}$, we might first normalize $\mathbf p^{\ast}$ and enter the normalized form into our sampling function. Here, we use an alternative approach and rather change the draw from the uniform distribution from the interval $(0,1)$ to $(0, s)$. As we are working on the log-scale, however, this means that we have to sum the exponentials of the log-probabilities and then take the logarithm again to obtain $s$ (which, again, can lead to underflow if some of the probabilities are small), because



$$ \log(s) = \log\left(\sum_{j=1}^K p_j\right)  = \log\left(\sum_{j=1}^K \exp[\log(p_j)]\right)$$



Fortunately, we have already calculate this quantity when accumulating the probabilities on the log-scale in the `lqvec` vector. That is, the last element of `lqvec` is `log(exp(lpvec[1]) + ... + exp(lpvec[K]))` which is exactly `log(sum(exp(lpvec)))`, the normalizing quantity we need.

Also, 



$$\Pr[-\log(sU)\le c] = 1 - \Pr[U \le  e^{-(c+\log s )}] = 1 - e^{-(c+ \log s)}, \quad c > 0$$



which is the distribution function of a $\text{shifted-Exponential}(1, -\log s)$ random variable. Thus, to sample uniformly from the interval $(0, s)$ on the log-scale, we might subtract from $\log s$ a draw from an $\text{Exponential}(1)$ distribution.

So, our final sampling function looks like the following:


{% highlight cpp %}
//[[Rcpp::export]]
arma::uword lsamp_one(const arma::vec &lpvec) {
    
    //check whether all elements are finite
    if (lpvec.has_inf())
        Rcpp::stop("log-probabilities have to be finite");

    // get indices of elements in descending order
    arma::uvec indx = arma::sort_index(lpvec, "descend");
    // logarithm of accumulated probabilities
    arma::vec lqvec = log_accu_exp(arma::vec(arma::sort(lpvec, "descend")));

    // sample log(unif(0,sum(exp(pvec))))
    double u = arma::conv_to<double>::from(lqvec.tail(1L));
    u -= R::rexp(1.0);

    // sample integer
    for (arma::uword k = 0; k < lqvec.n_elem; ++k) {
        if (u <= lqvec(k))
            return indx(k);

    }
    
    Rcpp::stop("could't find index (lsamp_one)");

}
{% endhighlight %}
Notice that this function will work for properly normalized probability vectors as well. In these situations, `lqvec.max()` would be equal to zero and so we have exactly the same sampling function as `lsamp_temp`.

We might also rewrite the function that accepts unnormalized probability vectors (instead of log-probability vectors) as follows:

{% highlight cpp %}
//[[Rcpp::export]]
arma::uword samp_one(const arma::vec &pvec) {
    
    // check that all elements are positive
    if (arma::any(pvec < 0.0))
      Rcpp::stop("probabilities have to be positive");

    // generate a vector of indices, so that 0 represents the largest 
    // element, 1 the secon largest, and so on
    arma::uvec indx = arma::sort_index(pvec, "descend");
    
    // generate the q-vector, the vector of cumulative sums
    arma::vec qvec = arma::cumsum(arma::sort(pvec, "descend"));
    
    // draw randomly from (0,s)
    double u = arma::conv_to<double>::from(qvec.tail(1L));
    u *= arma::randu();
    
    // find interval into which u falls
    for (arma::uword k = 0; k < pvec.n_elem; ++k) {
        
        if (u <= qvec(k))
            return indx(k);
        
    }
    
    Rcpp::stop("could't find index (samp_one)");
}
{% endhighlight %}


A small Monte Carlo simulation suggests that the code works well:

{% highlight r %}
set.seed(123)
# number of draws
n.sim = 5e5
# number of outcomes
k = 4

# generate (unnormalized) probability vector
p = runif(k)
lp = log(p)

# normalized probabilities
norm.p = p / sum(p)

# sample integers from unnormalized prob. vectors
res1 = replicate(n.sim, sample.int(k, 1L, prob = p) - 1L)
res2 = replicate(n.sim, samp_one(p))
res3 = replicate(n.sim, lsamp_one(lp))

# print results of simulation
print(
  rbind(
    norm.p, 
    samp.int  = prop.table(table(res1)),
    samp_one  = prop.table(table(res2)),
    lsamp_one = prop.table(table(res3))
  )
)
{% endhighlight %}

```
##                   0         1         2         3
## norm.p    0.1214495 0.3329164 0.1727188 0.3729152
## samp.int  0.1219100 0.3326480 0.1723480 0.3730940
## samp_one  0.1213660 0.3326760 0.1720160 0.3739420
## lsamp_one 0.1218680 0.3335720 0.1721560 0.3724040
```

The full `C++` code used is reproduced below:

{% highlight cpp %}
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]

inline double log_add_exp(const double x, const double y) {

    double d = x - y;
    if (d > 0.0)
        return x + std::log1p(std::exp(-d));
    if (d <= 0.0)
        return y + std::log1p(std::exp(d));

}

template <typename T>
inline T log_accu_exp (const T& x) {

    // initialize with x
    T y(x);
  
    // accumulate probabilities on the log-scale
    typename T::iterator i = y.begin() + 1;
    for (; i < y.end(); i++) {
        *i = log_add_exp(*(i - 1), *i);
    }

    return y;
}

//[[Rcpp::export]]
arma::uword samp_one(const arma::vec &pvec) {
    
    // check that all elements are positive
    if (arma::any(pvec < 0.0))
      Rcpp::stop("probabilities have to be positive");

    // generate a vector of indices, so that 0 represents the largest 
    // element, 1 the secon largest, and so on
    arma::uvec indx = arma::sort_index(pvec, "descend");
    
    // generate the q-vector, the vector of cumulative sums
    arma::vec qvec = arma::cumsum(arma::sort(pvec, "descend"));
    
    // draw randomly from (0,s)
    double u = arma::conv_to<double>::from(qvec.tail(1L));
    u *= arma::randu();
    
    // find interval into which u falls
    for (arma::uword k = 0; k < pvec.n_elem; ++k) {
        
        if (u <= qvec(k))
            return indx(k);
        
    }
    
    Rcpp::stop("couldn't find index (samp_one)");
    
}

//[[Rcpp::export]]
arma::uword lsamp_one(const arma::vec &lpvec) {
    
    //check whether all elements are finite
    if (lpvec.has_inf())
        Rcpp::stop("log-probabilities have to be finite");

    // get indices of elements in descending order
    arma::uvec indx = arma::sort_index(lpvec, "descend");
    // logarithm of accumulated probabilities
    arma::vec lqvec = log_accu_exp(arma::vec(arma::sort(lpvec, "descend")));

    // sample log(unif(0,sum(exp(pvec))))
    double u = arma::conv_to<double>::from(lqvec.tail(1L));
    u -= R::rexp(1.0);

    // sample integer
    for (arma::uword k = 0; k < lqvec.n_elem; ++k) {
        if (u <= lqvec(k))
            return indx(k);

    }
    
    Rcpp::stop("couldn't find index (lsamp_one)");

}
{% endhighlight %}