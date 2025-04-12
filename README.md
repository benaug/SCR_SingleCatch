# SCR_SingleCatch
Spatial capture recapture for single-catch traps using a latent ID approach. Also includes multi-catch traps observation model.

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

This approach is pretty slow. I did a simulation study with a N=50 (M=100), J=49, K=5, which I did not time
but I think it took about 1.5hr to run 20K iterations and could really get by with 10K. Run time will slow down raising
N, J, and K. I tried a few short runs with 50 x 144 x 5 and that wasn't too slow. Raising N does really slow things down. There
may be efficiency gains yet to be made. If not, perhaps this is a "smaller N" approach.

4/7/25: Added a version with a spatial density covariate. Realized N is random in the data generator, so data realizations with
larger realized Ns take longer to run.

4/8/25: Added a version with a behavioral response to capture. The assumption is that an individual has to be actually
captured, not "latently captured" to undergo a trap response. Therefore, the capture states (first/subsequent) are observed from
the observed capture history. Observation model parameters are estimated quite imprecisely, particularly first capture p0 and
sigma when the trap response is positive.

4/12/25: Added multi catch observation model from Efford and Borchers 2009 using a competing hazards model.

https://link.springer.com/chapter/10.1007/978-0-387-78151-8_11
Efford, Murray G., David L. Borchers, and Andrea E. Byrom. 
"Density estimation by spatially explicit capture–recapture: likelihood-based methods."
Modeling demographic processes in marked populations (2009): 255-269.

Added Mb + Dcov for single catch traps. Now, there is M0, Mb, M0 + Dcov, Mb + Dcov for single and multi catch traps. I have not tested the
multicatch versions. Will remove this disclaimer when tested.


