README

1. `main.m` contains routine for generating QoI's from prefactors
2. Line 24-26 of `main.m` are 3 example choices of the prefactors (theta)
3. to plot the PDE solutions, change Line 59 of `reactdiffuse1d2sp.m`



Data: MC simulation for 3 cases

* Case 1: `theta =[log(.01), log(0.4), 0.1,  -1,   0,   1, 0.9,   0,   0,  -1]`
  *  `w2.mat` likelihood err 20%, positive diffusivities NOT enforced, 30k iterations 
  * `w7.mat` likelihood err 15%, positive diffusivities, 30k iterations
  * `w8.mat` likelihood err 10%, 30k iterations, `w8cont` additional 30000 iters
* Case 2: `theta =[log(.03), log(1.2), 0.1,  -1,   0,   1, 0.9,   0,   0,  -1]`
  *  `w3.mat` likelihood err 20%, 30000 iterations, `w3cont.mat` additional 50000 iters (diverge!)
  * `w6.mat` likelihood err 10%, 30000 iterations, `w6cont.mat` additional 50000 iters, `w6cont2` another 30000 iters, `w6cont3` another 30000 iters (running)
* Case 3: `theta =[log( .1), log(  4), 0.1,  -1,   0,   1, 0.9,   0,   0,  -1]`
  * `w4.mat` likelihood err 20%, 30000 iterations (diverge!)
  * `w5.mat`  likelihood err 10%, 50000 iterations