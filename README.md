# SARS-CoV-2 NAAT and Antigen Testing Algorithms

These files can be used to reproduce the testing algorithms and simulations described in the article [**Quantitative Comparison of SARS-CoV-2 Nucleic Acid Amplification Test and Antigen Testing Algorithms: A Decision Analysis Simulation Model**](https://www.medrxiv.org/content/10.1101/2021.03.15.21253608v1). 

The R file [*AlgorithmFunctions.R*](https://github.com/CDCgov/SARS-CoV-2-NAAT-and-Antigen-Testing-Algorithms/blob/main/AlgorithmFunctions.R) contains the codes defining each algorithm described in the article:

<ol type="A">
  <li>NAAT Only</li>
  <li>Antigen Only</li>
  <li>NAAT Confirmation for Symptomatic Antigen-Negative  and Asymptomatic Antigen-Positive Results</li>
  <li>NAAT Confirmation of Negative Antigen Results </li>
  <li>Repeat Antigen Confirmation of Antigen-Negative Results</li>
  <li>NAAT for Asymptomatic Persons & Symptomatic Persons with Positive Antigen Results
</li>
</ol>

The R file [*AlgorithmSimulations.R*](https://github.com/CDCgov/SARS-CoV-2-NAAT-and-Antigen-Testing-Algorithms/blob/main/AlgorithmSimulation.R) calls *AlgorithmFunctions.R* to perform a specified number of simulations of each algorithm at different levels of prevalence.

The CSV file [*InputParameterRanges.csv*](https://github.com/CDCgov/SARS-CoV-2-NAAT-and-Antigen-Testing-Algorithms/blob/main/InputParameterRanges.csv) is used to specify the distribution function for each input parameter's random sample. It is also called by *AlgorithmSimulations.R*. 
