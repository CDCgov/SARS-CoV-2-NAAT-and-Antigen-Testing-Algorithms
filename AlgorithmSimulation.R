#Load code defining each algorithm
source("AlgorithmFunctions.R")

#Load/Install libraries for Latin Hypercube Sampling and Triangular distribution
if(!require("lhs")){install.packages("lhs")}
if(!require("triangle")){install.packages("triangle")}
library(lhs); library(triangle);

#Read CSV to parameterize triangular sampling distributions for each input parameter except prevalence (specified below)
param.in<-read.csv("InputParameterRanges.csv",stringsAsFactors = F)

#Specify the number of simulations to be performed for each algorithm
#Default is 50,000 simulations of each algorithm at each prevalence level
n.samples<-5e4

#Specify the prevalence levels to be evaluated 
#Default is to simulate each algorithm n.samples times at 4 prevalence levels: 5%, 10%, 15%, 20%
prevs<-seq(.05,.2,.05)

#Generate a Latin Hypercube Sample for n.samples simulations each with 12 input parameters (defined by rows of param.in)
#lhs.rand object will be a matrix object; each row corresponds to a unique set of 12 inputs randomly sampled across the range 0-1
lhs.rand<-randomLHS(n.samples,nrow(param.in))

#Prepare to transform LHS from uniform distribution ranging 0-1 for each parameter to triangular distribution with a specified point/range for each parameter
#Create a matrix to store transformed values
param.sample<-matrix(ncol=nrow(param.in),nrow=n.samples)
colnames(param.sample)<-param.in$param.name

#Loop over each parameter defined by param.in
for(i in 1:nrow(param.in)){
  #Use each randomly sampled value in lhs.rand (range 0-1) as a percentile value to draw from the triangular distribution defined by the parameterization in param.in
  #Use these transformed values from the triangular distribution to populate the corresponding column in param.sample
  param.sample[,paste0(param.in[i,"param.name"])]<-qtriangle(lhs.rand[,i],a=param.in[i,"low"],b=param.in[i,"high"],c=param.in[i,"point"])
  
}
#After looping, param.sample contains the input values for Monte Carlo simulation: each column corresponds to a different input parameter and each row corresponds to the input values for a unique simulation

#Create lists to store simulation results for each algorithm
#Each list will contain a data.frame for each prevalence level's simulation results
#Data.frames will have n.sample rows, each row containing the results of the algorithm corresponding to the same row of input values in param.sample

naat_only.res<-list() #Results from Algorithm (A) NAAT Only
ag_only.res<-list() #Results from Algorithm (B) Antigen (Ag) Only
naat_sx.agneg_asx.agpos.res<-list() #Results from Algorithm (C) NAAT Confirmation for Symptomatic Antigen-Negative (Sx/Ag-neg) and Asymptomatic Antigen-Positive (Asx/Ag-pos) Results
naat_agneg.res<-list() #Results from Algorithm (D) NAAT Confirmation of Negative Antigen Results (Ag-neg)
repeat_ag_agneg.res<-list() #Results from Algorithm (E) Repeat Antigen Confirmation of (Ag-neg)
naat_asx_sx.agpos.res<-list()#Results from Algorithm (F) NAAT for Asymptomatic Persons (Asx) & Symptomatic Persons with Positive Antigen Results (Sx/Ag-pos)

#Loop over each prevalence level to add data.frames to result lists
for(p in prevs){
  
  #Use apply() to iterate simulations across each row of param.sample
  #FUN= do.call() to call algorithm function, supplying parameters & prevalence as a list of inputs
  #apply() will return result of each simulation as a separate column; use t() to transform data.frame columns to rows (now each result row will correspond to a row of input values in param.sample)
  
  #Simulate Algorithm (A)
  naat_only.res[[paste0("p",p)]]<-t(apply(param.sample,1,function(x) do.call(naat_only.fun,as.list(c(prev=p,x)))))
  
  #Simulate Algorithm (B)
  ag_only.res[[paste0("p",p)]]<-t(apply(param.sample,1,function(x) do.call(ag_only.fun,as.list(c(prev=p,x)))))
  
  #Simulate Algorithm (C)
  naat_sx.agneg_asx.agpos.res[[paste0("p",p)]]<-t(apply(param.sample,1,function(x) do.call(naat_sx.agneg_asx.agpos.fun,as.list(c(prev=p,x)))))
  
  #Simulate Algorithm (D)
  naat_agneg.res[[paste0("p",p)]]<-t(apply(param.sample,1,function(x) do.call(naat_agneg.fun,as.list(c(prev=p,x)))))
  
  #Simulate Algorithm (E)
  repeat_ag_agneg.res[[paste0("p",p)]]<-t(apply(param.sample,1,function(x) do.call(repeat_ag_agneg.fun,as.list(c(prev=p,x)))))
  
  #Simulate Algorithm (F)
  naat_asx_sx.agpos.res[[paste0("p",p)]]<-t(apply(param.sample,1,function(x) do.call(naat_asx_sx.agpos.fun,as.list(c(prev=p,x)))))
  
  
}

