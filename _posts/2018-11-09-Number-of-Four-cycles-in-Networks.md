---
layout: post
title: Number of Four-cycles in Networks
date: 2017-06-26 21:48
author: baruuum
comments: false
categories: [Graph Theory, R]
---

I'm in the middle of moving my website from WordPress to Github Pages. I've exported my posts on WordPress into a `.xml` file and have used [this great program](https://github.com/theaob/wpXml2Jekyll) to transform the pages into Markdown documents. One problem that remains is translating equations from the blogs in Wordpress to those in Jekyll. So, I'm trying different things out. So, here is a test (by the way, the reason for working on this weird mathematical problem, if I remember correctly, was the amazing <em>AJS</em> article by Peter Bearman, James Moody, and Katherine Stovel on the sexual network of high school students).

<hr />

#### Number of Four-Cycles in a Graph

Let $ G=(V,E) $ be a graph with associated adjacency matrix $ A $. Consider a closed walk of length 4 starting at the vertex $ v_1\in V $. We have to consider four different scenarios:

- for each $ v_j $ adjacent to $ v_i $, there exists a length-four walk $ \{v_i, v_j,v_i,v_j,v_i\} $ which is obviously no 4-cycle. The number of such length-four walks are equal to the degree of $ v_i $;

- if $ \{v_i,v_j\}, \{v_j,v_k\} \in E(V) $, $ i\ne j\ne k $, then $ \{v_i,v_j,v_k,v_j,v_i\} $ is a walk of length 4 but no 4-cycle. The number of such walks is equal the degree of $ v_j $ minus one (the tie to $ v_i $);

- if $ v_i $ has more than one neighbor, then $ \{v_i,v_j,v_i,v_k,v_i\} $ is also a walk of length 4, but not a cycle. For a vertex with more than one neighbor, the number of these walks are given as $ 2\binom{deg(v_i)}{2} $;

- if $ \{v_i,v_j,v_k,v_l,v_i\} $ is a length-4 closed walk visiting distinct vertices, then so is $ \{v_i,v_l,v_k,v_j,v_i\} $. Hence we are always double counting.


Let $ N=\{i: v_i \in V(G)\} $ be an index set and denote the degree of vertex $ v_i $ as $ d_G(v_i) $. If we let $ c_i $ be the number of 4-cycles in which node $ i $ is involved, then

$$  \begin{aligned} a_{ii}^{(4)} &= 2 c_i + d_G(v_i) + \sum_{j\in N}a_{ij}\Big(d_G(v_j)-1\Big) + 2\binom{d_G(v_i)}{2} \\ &= 2 c_i + d_G(v_i) + \sum_{j\in N}a_{ij}\left(\sum_{k\in N}a_{jk} -1\right) + 2\binom{d_G(v_i)}{2} \\ &= 2 c_i + \sum_{k\in N}a_{ik}^{(2)} + 2\binom{d_G(v_i)}{2} \end{aligned}  $$

as $ d_G(v_i) = \sum_{j\in N}a_{ij} $. Hence,

$$  \begin{aligned} c_i&= \frac{1}{2}\left(a_{ii}^{(4)} - \sum_{j\in N}a_{ij}^{(2)} - 2\binom{d_{G}(v_i)}{2}\right). \end{aligned}  $$

Since each cycle of length 4 consists of four vertices, it follows that the number of 4-cycles $ |\mathcal C_4| $ in $ G $ is

$$  \begin{aligned} |\mathcal C_4(G)| &= \frac{1}{4}\left[\sum_{i\in N} \frac{1}{2}\left(a_{ii}^{(4)} - \sum_{j\in N}a_{ij}^{(2)} - 2\binom{d_{G}(v_i)}{2}\right)\right] \\ &= \frac{1}{8}\left[\sum_{i\in N} a_{ii}^{(4)} - \sum_{i,j\in N}a_{ij}^{(2)} - 2\sum_{i\in N}\binom{d_{G}(v_i)}{2}\right]. \end{aligned}  $$

Now, we notice that the first term is the trace of the matrix $ A^4 $ and that $ \sum_{i\in N} d_G(v_i) = 2|E(G)| $. Further,

$$  \begin{aligned} \sum_{i,j\in N}a_{ij}^{(2)} &= \sum_{i\in N} a_{ii}^{(2)}+ 2\sum_{i\ne j \land i,j \in N}a_{ij}^{(2)} \\ &=\sum_{i\in N} d_G(v_i) + 2\sum_{i\in N}\binom{d_G(v_i)}{2}, \end{aligned}  $$

i.e., the sum of the elements in $ A^2 $ is equal to the diagonals plus the off-diagonals. But the sum of the diagonals is equal to the sum of the degree sequence of the graph and the sum of the off-diagonals is equal to the number of length-2 paths. The number of length-2 paths can be calculated, in turn, as follows: choose a vertex, say $ v_i $; if the vertex has a degree greater than one, calculate the number of ways in which any two neighbors of $ v_i $ can be chosen, each of which determines a unique length-2 path for which $ v_i $ lies at the center. Summing over this number for all vertices in the graph gives the number of length-2 paths, as no two different vertices can lie at the center of the same length-2 path. On the other hand, if $ d_G(v_i) \le 1 $, $ \binom{d_G(v_i)}{2}=0 $.

Therefore it follows that

$$  \begin{aligned} |\mathcal C_4(G)|&= \frac{1}{8}\left( \text{trace}(A^4) - 2|E(G)| - 4\sum_{i}\binom{d_G(v_i)}{2}\right). \end{aligned}  $$

<hr />

Works great! The `R` code that was used is shown below:

{% highlight r %}
library("magrittr")

# file names
fn <- "filname.md"
out <- "output.md"

# convert equations
x <- readLines(fn) %>%
    gsub("&amp;fg=0*", "", .) %>%
    gsub("&amp;=&amp;", "&=", .) %>%
    gsub("(<\\s*p.*?>)(.*?)(<\\s*/p\\s*>)", "\\2", .) %>%
    gsub("(\\$\\s*latex\\s*\\\\displaystyle)(.*?)(\\$)", 
         "\n\n$$ \\2 \n\n", .) %>%
    gsub("\\\\begin\\{array\\}\\{.*?\\}", "\\\\begin{aligned}", .) %>%
    gsub("\\\\end\\{array\\}", "\\\\end{aligned}", .) %>%
    gsub("(\\$\\s*latex\\s*\\{)(.*?)(\\}\\s*\\$)", "$$ \\2 $$", .)

# check
cat(x)

# write out
writeLines(x, out)
{% endhighlight %}

