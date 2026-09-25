# wmdc

**Weighted Marginal Distribution Calibration**

An R package for calibrating impedance functions used in spatial analysis. 
Given the marginal distribution of an impedance variable (e.g. travel time or 
distance) expressed as a frequency table or vector of observations, `wmdc` 
estimates the parameters of a chosen impedance decay function using maximum 
likelihood estimation. The package currently supports exponential, Gaussian, and 
power-exponential forms, and includes a built-in spatial weight function to
correct for the systematic variation in the structure of origins and 
destinations over impedance, which is a source of bias and misspecification in 
many existing implementations.

## Installation

```r
# install.packages("pak")
pak::pak("harrysroberts/wmdc")
```

## Background

In spatial analysis, an **impedance function** describes how willingness to 
travel to a destination decays with the increasing cost of reaching it, 
expressed in terms of an **impedance** variable (distance, travel time, etc.). 
This is a key component of both spatial interaction models, used to predict 
interzonal flows, and accessibility measures, in which the impedance function is 
used to weight the contribution of each destination to the overall accessibility 
of an origin location.

A common method of calibrating impedance functions, particularly in 
accessibility analysis, is fitting to the observed marginal distribution of the 
impedance variable. For example, if travel time is taken to be the impedance 
variable, then the distribution of travel times in the study area can be used to 
estimate the parameters of an impedance function.

However, most existing implementations fail to account for  the underlying
spatial structure constraining travel behaviour, which can lead to biased 
estimates of the decay parameters and misspecification of the impedance function. 
For instance, the relative lack of short-distance trips in observed data is often
interpreted as evidence of a 'frictionless' behavioural response to impedance at 
short distances, thus supporting a Gaussian-type decay function 
(e.g. [Ingram 1971](https://doi.org/10.1080/09595237100185131)). A more likely 
explanation for this pattern is that there are simply fewer reachable 
destinations at short distances, which is a property of the spatial structure 
rather than a behavioural response to impedance.

To address this issue, `wmdc` models the observed marginal distribution of the
impedance variable $p(x)$ as a product of two components: 

- the **impedance function** $f(x)$ itself, capturing the behavioural response 
to impedance; and

- a **spatial weight function** $g(x)$ that captures how the effect of spatial
structure varies with impedance, and is closely related to the number of 
reachable destinations at impedance $x$.

Thus the marginal impedance distribution is expressed as the product 
$p(x) \propto g(x) \cdot f(x)$.

In this package, the impedance function is modelled using a flexible family of
power-exponential functions

$$f(x) = e^{-(\lambda x)^r}$$

where $\lambda$ is a rate parameter and $r$ controls the shape of the decay.
Exponential decay corresponds to $r = 1$, Gaussian decay to $r = 2$, and 
power-exponential decay to the general case, where $r$ is estimated from the 
empirical data.

The spatial weight function is modelled as a power function of the impedance 
variable

$$g(x) \propto x^{\alpha - 1}$$

for some $\alpha \geq 1$. In the case $\alpha = 2$, this is linear, reflecting the 
fact that if destinations are distributed uniformly in a 2-D plane, the number 
of destinations at distance $x$ for a given origin is proportional to the 
circumference of the circle of radius $x$ around that origin. In $\alpha = 1$
and the count represents an uncorrected model with constant weighting, and 
corresponds to a 1-D spatial structure in which the number of destinations at 
distance $x$ is constant. In general $\alpha$ is related to the effective 
spatial dimension of the study area, which may be non-integer due to the fractal
nature of transport networks.

Taking the product of the impedance and spatial weight functions and normalising
gives the probability density function for a **generalised gamma distribution** 

$$p(x) = \frac{r \lambda^\alpha}{\Gamma(\alpha/r)} 
x^{\alpha-1} e^{-(\lambda x)^r}$$

where $\Gamma(\cdot)$ is the gamma function. The free parameter(s) of this
distribution are estimated using maximum likelihood estimation over the data.

The impedance function with the calibrated $\lambda$ and $r$ parameters
(although $r$ may be fixed) can then be used in spatial interaction and 
accessibility models. The spatial weight function, density function, and 
survival function are also provided as outputs, but are not intended for use in
spatial interaction models. In particular, the output density and survival 
functions **should not be misinterpreted as the impedance function itself** as 
they are inclusive of the spatial weight function, which captures the spatial 
structure rather than travel behaviour.

## Usage

The main function is `wmdc()`. It takes a two-column data frame of 
impedance values and their frequencies, and returns a list of calibrated 
functions and diagnostics.

```r
library(wmdc)

# Example: travel time frequency table
travel_times_table <- data.frame(
  time    = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60),
  count   = c(100, 200, 250, 200, 100, 200, 60, 60, 70, 40, 20, 60)
)

result_1 <- wmdc(travel_times_table)

# Example: vector of travel time observations
travel_times_vector <- c(
    5.7, 7.3, 7.5, 8.9, 9.7, 10.9, 11.1, 13.7, 14.5, 18.4, 
    21.5, 24.1, 27.5, 27.7, 28.9, 32.8, 38.3, 38.5, 43.9, 52.0
    )

result_2 <- wmdc(travel_times_vector)

```

### Arguments

| Argument | Default | Description |
|---|---|---|
| `data` | — | A two-column data frame: impedance values and frequencies |
| `shape` | `1` | Fixed shape of impedance function ($r$) (or starting value if `free_shape = TRUE`)|
| `dimension` | `2` | Fixed effective spatial dimension ($\alpha$) (or starting value if `free_dimension = TRUE`)|
| `free_shape` | `FALSE` | If `TRUE`, estimate shape $r$ as a free parameter|
| `free_dimension` | `FALSE` | If `TRUE`, estimate dimension $\alpha$ as a free parameter|
| `maxLik_output` | `FALSE` | If `TRUE`, include full `maxLik` optimisation output in results|

### Impedance function forms

| Form | `shape`| `free_shape` | Impedance function |
|---|---|---|---|
| Exponential | 1 | FALSE | $f(x;\lambda) = e^{-\lambda x}$ |
| Gaussian | 2 | FALSE | $f(x;\lambda) = e^{-(\lambda x)^2}$ |
| Power exponential (fixed) | Fixed $c>0$| FALSE | $f(x;\lambda) = e^{-(\lambda x)^c}$ |
| Power exponential (free) | Starting value | TRUE | $f(x;r,\lambda) = e^{-(\lambda x)^r}$ |


### Return value

`wmdc()` returns a named list:

| Element | Description |
|---|---|
| `impedance_function` | Calibrated impedance decay function $f(x)$ |
| `spatial_weight_function` | Calibrated spatial weight function $g(x)$ |
| `density_function` | Probability density function $p(x)$ for the marginal impedance distribution |
| `survival_function` | Survival function $P(X \geq x)$ |
| `parameters` | Named vector of calibrated parameter estimates for impedance function |
| `parameters_se` | Standard errors of parameter estimates for impedance function |
| `dimension` | Effective spatial dimension $\alpha$ (fixed or estimated) in spatial weight function |
| `dimension_se` | Standard error of estimated dimension (`NA` if `free_dimension = FALSE`) |
| `log_likelihood` | Log-likelihood at the maximum |
| `AIC` | Akaike Information Criterion |
| `survival_R_sq` | $R^2$ for the survival function fit |
| `maxLik_output` | Full `maxLik` output (if `maxLik_output = TRUE`) |

### Example

```r
result <- wmdc(travel_times)

# Inspect parameters
result$parameters
result$parameters_se

# Model fit
result$AIC
result$log_likelihood
result$survival_R_sq

# Apply the calibrated impedance function
x <- seq(0, 60, by = 1)
fx <- result$impedance_function(x)
plot(x, fx, type = "l", xlab = "Travel time (min)", ylab = "Impedance")
```

## FAQs

**Is WMDC accurate?**

The weighting in the WMDC method is an approximation to the systematic variation 
in the spatial structure over impedance, capturing growth in the number of 
reachable destinations with increasing impedance in an $\alpha$-dimensional 
space. This may be an improvement over the unweighted method, but it does not 
account for all sources of distortion arising from the spatial structure.

In particular, it struggles in situations where the study area is small and the 
impedance decay is weak, as interference from the boundary dominates the 
marginal distribution and leads to systematic overestimation of the decay rate.

Moreover, where the decay is strong but the study area is sparsely populated, 
trip patterns are more deterministic, hindering accurate calibration (although 
scenarios like these are unrealistic and rare in the real world).

Thus, WMDC is most accurate when the study area is large relative to the half 
life of impedance decay (e.g. over a national or regional scale), provided that
there is sufficient density.

For more detail, see [`here`](https://github.com/harrysroberts/hEART_2026).

**What does the dimension parameter $\alpha$ mean, and what values should it take?**

As noted above, $\alpha$ is related to the effective spatial dimension of the 
study environment. This is a fractal extension of the concept of spatial
dimension, which can take non-integer values. Given that transport networks span
a subset of the two-dimensional plane, it is expected that $\alpha$ should fall 
between 1 and 2.

Therefore, if setting `free_dimension = TRUE`, it is advisable to check that the 
estimated $\alpha$ is within this range. If the estimated $\alpha$ falls outside 
this range, the calibration should be re-run with `free_dimension = FALSE` and a 
fixed value of `dimension` between 1 and 2 inclusive.

**Will the calibration always produce the optimal solution?**

If either the shape parameter $r$ or the dimension parameter $\alpha$ is 
estimated as a free parameter, through selecting `free_shape = TRUE` or
`free_dimension = TRUE` respectively, there is no guarantee of a unique solution
to the maximum likelihood estimation problem, as the Hessian matrix of the 
log-likelihood function may not be negative definite. If the optimisation fails 
to converge, or converges to a solution with parameters outside the expected 
parameter space, it is advisable to re-run the calibration with fixed 
parameters, or to try different starting values for the optimisation.


**Why is an exponential impedance function used as the default?**

The exponential form of the impedance function is used as the default in this
package (through `shape = 1`). This is not an arbitrary decision: exponential 
decay emerges as the entropy-maximising solution to trip distribution problem 
under the minimal assumption of an additive impedance budget. Put simply, if all
one assumes about the behaviour of travellers is that they are are limited in 
the amount of distance or time they are willing to travel, then an exponential 
impedance function is the 'best guess'. Using this as the default is therefore 
a principled choice, consistent with Ockham's razor.

Users may wish to use a different form of the impedance function, such as 
Gaussian or power-exponential, for two main reasons:

- They have a different hypothesis about the behavioural response to impedance
- They wish to achieve a better fit to the observed data

In the latter case, users are advised to be mindful of whether the improvement
over the exponential form is substantive, to avoid overfitting the data and 
losing the parsimony of the exponential form. Additionally, they are encouraged
to consider the behavioural implications of the shape of the impedance function,
as deviation from the exponential form suggests alternative or additional 
assumptions about travel behaviour. While 'letting the data speak for itself' 
has become a popular paradigm, reflection on these implications may elicit 
valuable insights, or bring to light the existence of structural bias in the 
data that may be causing it to misspeak.

Users wishing to learn more about the relationship between impedance functions
and entropy are encouraged to read [Wilson (1967)](https://doi.org/10.1016/0041-1647(67)90035-4) 
and [Anas (1983)](https://doi.org/10.1016/0191-2615(83)90023-1).

## Dependencies

- [`maxLik`](https://cran.r-project.org/package=maxLik) — maximum likelihood estimation

## Citation

If you use this package in your research, the following citations are appreciated:

**Poster presented at hEART 2026:**

 Roberts, H. S., Calastri, C., Batley, R. 2026. _Using marginal impedance distributions to calibrate 
impedance functions for accessibility measurement_ \[Poster\]. European Association for Research in Transportation (hEART), 2026, Paris.

**This software:**

Roberts, H.S. (2026) “wmdc”. Zenodo. doi:TBC.

## License

MIT License. See LICENSE file for details.

## Author

Harry Roberts ([H.S.Roberts@leeds.ac.uk](mailto:H.S.Roberts@leeds.ac.uk))

Institute for Transport Studies, University of Leeds

