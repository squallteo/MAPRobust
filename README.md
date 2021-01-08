This repository contains the R programs accompanying the manuscript: 

Zhang, H., Chiang, A., Branson, M., **On the Implementation of Robust Meta-Analytical-Predictive Prior**, *Statistics in Biopharmaceutical Research*

**Instructions**

1. Files starting with "01" are for the normal endpoint, whereas those starting with "02" are for the binary endpoint.  

2. The files "XXXSimRun.R" run the simulations. 

3. "XXXSimSpec0.R" contains the historical data and simulation specs corresponding to the null scenario (no treatment effect), whereas "XXXSimSpec1.R" to the alternative scenario. 

4. The simulation stores a series of R workspace files. The can be summarized by "XXXResultsSummary.R" files. 

5. "XXXDataAnalysis.R" files perform data analysis of one set of simulated current trial data. It also generates the plots in the manuscript. 

6. The .bugs files are model specifications in BUGS language. They must be stored locally and called by the file path. **They cannot be sourced in R** because they specify the half-normal priors using syntax that R doesn't recognize. 