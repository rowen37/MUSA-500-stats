# HW4 Part B report draft (Ripley's K-function)

Numbers below come from `code/HW4_Kfunction_analysis.R` → `output/table_02_kfunction_bands.csv`
and `output/fig_04_kfunction_envelope.png`. Sourced from the shared
`data/processed_prepare/markets_ppp.Rds` and `city_owin.Rds`, seed 671, 99 simulations,
identical to Part A.

## Methods (K-function section, drop-in text)

Ripley's K function asks a cumulative question: given one farmers market, how many other
markets fall within distance *d* of it, compared with what a random pattern would give?
Formally, K(d) = (1/λ)·E[number of other events within d of an arbitrary event], where λ is
the average intensity, here 62 markets over 3,970,679,605 square feet, or 1.561e-08 per
square foot. Under complete spatial randomness (CSR), which assumes a constant intensity and
independent locations, K(d) equals πd². Because the raw K curve grows as a parabola even
under randomness, we also report the variance-stabilised transform
L(d) = √(K(d)/π) − d, where CSR sits exactly at zero and positive values mean clustering.

We estimate K with the unbiased edge-corrected estimator
K̂(d) = (A/n²)·ΣᵢΣⱼ≠ᵢ I(dᵢⱼ ≤ d)·wᵢⱼ, where A is the study-area area, n the number of
markets, and wᵢⱼ is Ripley's isotropic edge weight that inflates pairs near the city
boundary, where part of the circle of radius dᵢⱼ falls outside the window and would
otherwise go uncounted. This is the R counterpart of ArcGIS's boundary correction. We did
not use *Simulate Outer Boundary Values* because the simulation-based correction requires
an outer buffer around the study polygon and a second round of simulated patterns outside
it; estimating K directly with the isotropic weight keeps the observed curve and all 99
reference curves on exactly the same footing inside one window, which is what the envelope
comparison needs.

Distance bands follow the assignment settings: beginning distance 0, increment 2,500 feet,
ten bands, so the outer edge is 25,000 feet (4.73 miles). The increment is defensible for
this city: 2,500 feet is about half a mile, roughly the walk or short drive a shopper would
accept, and ten bands reach far enough to cover the whole built-up extent of Philadelphia
without pushing the edge correction past its useful range.

Inference is by Monte Carlo. We generated 99 CSR patterns with the same number of points in
the same Philadelphia window (seed 671), re-estimated K for each with the same correction,
and built a pointwise envelope from the minimum and maximum of the 99 curves. A distance is
flagged when the observed K leaves that envelope; the corresponding one-sided p-value is the
rank of the observed value among the 99 simulated values, (1 + #extreme)/(99 + 1), so the
smallest reportable p is 0.010. To avoid reading ten separate tests as ten independent
chances to find something, we also ran a global test: for each pattern we took the largest
standardised departure of K from the simulated mean across all distances, and compared the
observed maximum with the distribution of that maximum over the 99 simulations.

A nonhomogeneous version of K replaces the constant λ with an intensity function λ(u)
estimated from an external reference measure, typically population, and compares the
observed pattern with an inhomogeneous Poisson process rather than a uniform one. That
version is needed when the null should be "markets follow people" instead of "markets are
spread evenly over the city polygon", which matters here because the Philadelphia boundary
includes the airport, the port, Fairmount Park and other tracts where nobody lives. We did
not run it: the assignment supplies zip-level median income, not a population surface fine
enough to define λ(u) at the market scale, and a zip-level intensity would smear the very
within-city variation the test is meant to detect.

## Results (K-function section, drop-in text)

The observed K curve lies above the 99-simulation envelope at all ten band edges, from
2,500 feet out to 25,000 feet (table_02, fig_04). At 2,500 feet K is 37.8 million square
feet against an envelope of 8.4 to 31.8 million; at 25,000 feet it is 3.11 billion against
1.78 to 2.26 billion. The L transform makes the same point in interpretable units: L(d) is
positive at every band, rises to a peak of about 7,260 feet at d = 20,000 feet (3.79 miles),
then eases to 6,438 feet at the outer edge. Read plainly, a market has on the order of
7,000 feet, or about 1.4 miles, more neighbouring markets than a random pattern would put there,
and the excess is present at every scale we tested, strongest around the three-to-four mile
range.

The pointwise rank p-value is 0.010 at every band edge from 2,500 feet outward, the smallest
value 99 simulations can produce. The global test agrees: the largest standardised departure
is 18.4, and none of the 99 simulated patterns came close, giving p = 0.010. So we reject
CSR for the Philadelphia farmers market pattern, and the direction of the departure is
clustered, not dispersed.

This matches Part A. The nearest-neighbour analysis gave R = 0.778 with the course formula
and R = 0.709 with the CDF-corrected Clark-Evans test, both well below 1, both rejecting CSR
at α = 0.05, and the Monte Carlo version returned p = 0.010. Two methods with different
windows of sensitivity — nearest neighbour looks only at the closest pair, K accumulates
pairs at every distance up to 25,000 feet — land on the same conclusion.

Two caveats belong with these numbers. First, the edge correction matters at this scale:
without it, K at 25,000 feet falls to 2.30 billion instead of 3.11 billion, roughly a
quarter lower, because pairs whose circles cross the city boundary are undercounted
(fig_05). The corrected curve is the one to interpret. Second, the 25,000-foot band sits
just beyond spatstat's recommended maximum of about 22,400 feet for this window; the
isotropic weight is still defined there and the curve is still monotone, but the outermost
band should be read as indicative rather than precise. Finally, the coincident pair at
Reading Terminal Market (Reading Terminal Market and Fair Food Farmstand @ RTM share one
address) is kept in every run, as agreed with Part A; at d = 0 it contributes a single
zero-distance pair, which is visible as the small non-zero K at the origin in the
uncorrected curve and does not change the conclusion at any band edge.
