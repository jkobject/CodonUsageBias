# Remarks about the code and what I did there


## utils

- because of hyper optimised C function it might is faster to run some unoptimized function relying heavily on C libraries (pandas, numpy, scipy...) than to try to optimize one's algorithm in python (needs to be tested with %%timeit)

- the multivariate normal approximation to the multinomial distribution does not seem to work if I provide a set of probabilities that are the same for each draw dimension because the matrix is not invertible since it is not SDP. SO the function might only work given non-symetrical probailities

- computejerem is faster than any other partition function and stays faster than even the randomdraw until reaching approx 30 000 000 possibilities

- trying to compute the complexity of the regular algorithm with for loops of for loops, I did not managed to do it for nbcod = 6 (see paper document on the subject)

- computepartition_sorted_full

- the random draw uses heavily dynamic programming and coding tricks to speed itself up. one of the trick is to 

- overall, the compute partition functions have allowed me to explore HPC and programmation issue when doing such highly repetitive computation on python as well as the many tricks to improve performance and speed. A great exploration of computational complexity

-getloc is really about recoding what someone wrote and enabled me to really fully understand the code of Yun Deng, it uses also heavily dynamic programming


- because of the way I do my save, basically storing everything into a dict and then a json, I have quite a few problem (typing issues, memory needed (approx 5times more than the final output file..))


#statistics about the code

- the code is about 6500 lines of pure code 
- it is 90% python, 5% javascript, 3% bash, 2% R and it contains 13 additional python packages and 2 additional R packages
 