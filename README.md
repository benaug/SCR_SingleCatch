# SCR_SingleCatch
Spatial capture recapture for single-catch traps using a latent ID approach

This is a latent ID SCR approach to the single catch trap model of Efford 2004 Oikos. Efford uses inverse prediction for estimation, 
whereas this approach is Bayesian.

https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/j.0030-1299.2004.13043.x

Efford, Murray. "Density estimation in live‐trapping studies." Oikos 106.3 (2004): 598-610.

See simulation algorithm on the final page.

Also, see Efford 2023 MEE

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14088

Efford, Murray G. "ipsecr: An R package for awkward spatial capture–recapture data." Methods in Ecology and Evolution 14.5 (2023): 1182-1189.

This Bayesian approach models the "true" capture history that would have occurred if traps did not fill up.
To do so, we need a "thinning model" to explain the capture history we observe conditional on the true capture
history. The thinning distribution here is exactly the same as the simulation algorithm. Capture times in this 
model are assumed to be exponential and capture events occur when an exponential RV is less than 1. Conditioning
on the true capture history on each occasion, y.true[1:M,1:J,k], the capture times of captured individuals are
exponential RVs right-truncated at 1. Then, we can use numerical integration to approximate the probability 
one exponential RV right-truncated at 1 is less than one or more exponential RVs right-truncated at one.
Compare dThin to the data simulator to explore further.

The MCMC approach updates the latent capture history y.true, and the latent capture order for the observed capture events.
For trap-occasions with observed captures, we need to update y.true[1:M,j,k] for other individuals. For trap-occasions with no
observed capture, we need to update the y.true[i,1:J,k] to consider the possibility they would have been captured at other
traps on this occasion if they were not stuck in the trap where they were actually captured. Finally, the capture order for
the observed captures must be updated. This is all in "ySampler" in the NimbleFunctions file.

This approach is pretty slow, but faster than I expected. I did a simulation study with a M=100, J=49, K=5, which I did not time
but I think it took about 1.5hr to run 20K iterations and could really get by with 10K. Run time will slow down raising
M, J, and K. I tried a few short runs with 100 x 144 x 5 and that wasn't too slow.

At some point, I'll add Mb.