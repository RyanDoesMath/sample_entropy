# sample_entropy
This repo houses the code for a paper we are writing that uses sample entropy on an open source dataset called VitalDB.  
The reason this code exists is because the python implementation on the Wikipedia page for sample entropy is (1) about 300 times slower than this code, (2) can't **really** be parallelized in python, and (3) is wrong.  

## What is Sample Entropy?  
Sample entropy (often abbreviated as sampen) is a metric for quantifying the uncertainty or stochasticity of time series data. If data varies in new an previously unseen ways, it has a higher value of sample entropy, and if the data repeats itself often, it will have a low value of sample entropy. It is derrived from approximate entropy, but removes a bias from the computation that approximate entropy suffers from (namely, self counting when determining matches). 
