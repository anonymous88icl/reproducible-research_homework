# Reproducible research: version control and R


The fourth and fifth questions are answered in this markdown file
developed using Quarto.

------------------------------------------------------------------------

## Questions 1-3

Questions 1 to 3 are answered in this linked [GitHub
Repository](https://github.com/anonymous88icl/logistic_growth/).

------------------------------------------------------------------------

## Question 4

### a)

#### Figure 1: Result of Executing `random_walk.R`

![](README_files/figure-commonmark/unnamed-chunk-1-1.png)

As we can see, the code produces two graphs that are the result of
random walks. In each of the 500 time-steps, the program generates an
uniformly random angle between $0$ and $2\pi$ and walks in that
direction by 0.25 units. The program records the coordinates at each
time step and eventually outputs the path in a graph. The time is
indicated by the color of the path (gets lighter as time passes). One
interesting thing we observe here is that the two graphs are different
even though they are produced by the same function. Also interestingly,
if we re-run the program, we will also get different results. This is
not ideal if we want reproducible results. (I didn’t make any specific
observations about the paths because that is also generated randomly.)

### b)

A random seed is an initial value used to initialize a random number
generator. Since a computer program is always fixed, it is impossible
for it to generate numbers in a “perfectly random” fashion. You could
always use something like time in a nanosecond-scale, but even then you
can’t guarantee complete randomness. Hence, the general principle to
generate these pseudorandom numbers is to start with an initial seed,
which we can either fix or set as some function of time (to make it more
random). Then, we apply some complex mathematical process with no simple
pattern to generate a sequence of numbers that seam to be random.

Hence, if we fix the seed, the outputs will still be close to random,
what it will be the same every time we re-initialize the random
function. We will now add this change to our code to make our random
walk reproducible.

### c)

We set the random seed to be equal to 880088, and the following output
is produced:

![](README_files/figure-commonmark/unnamed-chunk-2-1.png)

------------------------------------------------------------------------

## Question 5

------------------------------------------------------------------------
