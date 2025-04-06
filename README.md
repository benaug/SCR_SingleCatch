# SCR_SingleCatch
Spatial capture recapture for single-catch traps using a latent ID approach

This is a latent ID SCR approach to the single catch trap model of Efford 2004 Oikos. Efford uses inverse prediction for estimation, 
whereas this approach is Bayesian.

https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/j.0030-1299.2004.13043.x
Efford, Murray. "Density estimation in live‐trapping studies." Oikos 106.3 (2004): 598-610.

See simulation algorithm on the final page.

Also, see Efford 2023 MEE

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14088

This Bayesian approach models the "true" capture history that would have occurred if traps did not fill up.
To do so, we need a "thinning model" to explain the capture history we observe conditional on the true capture
history. The thinning distribution here is exactly the same as the simulation algorithm. Capture times in this 
model are assumed to be exponential at capture events occur when an exponential RV is less than 1. Conditioning
on the true capture history on each occasion, y.true[1:M,1:J,k], the capture times of captured individuals are
exponential RVs right truncated at 1. Then, we can use numerical integration to approximate the probability 
one exponential RV right truncated at 1 is less than one or more exponential RVs right truncated at one.
Compare dThin to the data simulator to explore further.

This approach is pretty slow. I did a simulation study with a M=100, J=49, K=5, which wasn't too slow, but 
I expect raising M, J, and K will slow things down pretty quickly.