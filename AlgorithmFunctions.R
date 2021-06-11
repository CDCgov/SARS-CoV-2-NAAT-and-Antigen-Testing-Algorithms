#Define Input Parameters

  #sx.cases : Proportion of infected population reporting symptoms at the time of testing 
  #sx.noncases : Proportion of uninfected population reporting symptoms at the time of testing
  #ag.sens.sx : Sensitivity of antigen test in symptomatic infected population
  #ag.sens.asx : Sensitivity of antigen test in asymptomatic infected population
  #ag.spec.sx  : Specificity of antigen test in symptomatic uninfected population
  #ag.spec.asx : Specificity of antigen test in asymptomatic uninfected population
  #naat.sens : Sensitivity of NAAT
  #naat.spec  : Specificity of NAAT
  #ag.sens.repeat : Sensitivity of a repeated antigen test in the infected population, conditional on an initial negative antigen result
  #ag.spec.repeat : Specificity of a repeated antigen test in the uninfected population, conditional on an initial negative antigen result
  #prev : Prevalence of infection in the population
  #asx.noncases.contacts : Proportion of asymptomatic uninfected population reporting having a recent close contact with infected person
  #naat.delay : mean days elapsed between sample collection and return of NAAT result
  

#Define Algorithm Model Outputs

  #detected_cases : Number of Infected persons correctly identified per 100,000 seeking evaluation
  #missed_cases : Number of Infected persons incorrectly identified (determined to be "uninfected" at the end of evaluation) per 100,000 seeking evaluation
  #false_positives : Number of Uninfected persons incorrectly identified (determined to be "infected" at the end of evaluation) per 100,000 seeking evaluation
  #n_ag_tests : Number of antigen tests performed per 100,000 seeking testing
  #n_naat_tests : Number of NAATs performed per 100,000 seeking testing
  #days_lost_waiting : Number of person-days of lost productivity among uninfected persons assumed to quarantine while awaiting results of NAATs


#Define Algorithms
  
  #(A) NAAT Only - each person is tested for SARS-CoV-2 infection by a NAAT
  naat_only.fun<-function(
    sx.cases = NA,
    sx.noncases  = NA,
    ag.sens.sx = NA,
    ag.sens.asx = NA,
    ag.spec.sx = NA,
    ag.spec.asx = NA,
    naat.sens = NA,
    naat.spec  = NA,
    ag.sens.repeat = NA,
    ag.spec.repeat = NA,
    prev = NA,
    asx.noncases.contacts = NA,
    naat.delay = NA
    
  ){
    detected_cases = 100000* prev*naat.sens
    missed_cases = 100000* prev*(1- naat.sens)
    false_positives= 100000*(1- prev)*(1- naat.spec )
    n_ag_tests = 0
    n_naat_tests = 100000
    days_lost_waiting = 100000*(1- prev)*(sx.noncases +  (1- sx.noncases)*asx.noncases.contacts )*naat.delay
    
    out=c(detected_cases,missed_cases,false_positives,n_ag_tests,n_naat_tests,days_lost_waiting)
    names(out)<-c("detected_cases","missed_cases","false_positives","n_ag_tests","n_naat_tests","days_lost_waiting")
    return(out)
    
  }
  
  #(B) Antigen (Ag) Only - each person is tested using a single antigen test, the result of which is used as a definitive diagnosis.
  ag_only.fun<-function(
    sx.cases = NA,
    sx.noncases  = NA,
    ag.sens.sx = NA,
    ag.sens.asx = NA,
    ag.spec.sx = NA,
    ag.spec.asx = NA,
    naat.sens = NA,
    naat.spec  = NA,
    ag.sens.repeat = NA,
    ag.spec.repeat = NA,
    prev = NA,
    asx.noncases.contacts = NA,
    naat.delay = NA
    
  ){
    detected_cases = 100000*prev*( sx.cases*ag.sens.sx + (1-sx.cases)*ag.sens.asx )
    missed_cases = 100000*prev*( sx.cases*(1-ag.sens.sx) + (1-sx.cases)*(1-ag.sens.asx) )
    false_positives= 100000*(1-prev)*( sx.noncases*(1-ag.spec.sx) + (1-sx.noncases)*(1-ag.spec.asx) )
    n_ag_tests = 100000 
    n_naat_tests = 0
    days_lost_waiting = 0
    
    out=c(detected_cases,missed_cases,false_positives,n_ag_tests,n_naat_tests,days_lost_waiting)
    names(out)<-c("detected_cases","missed_cases","false_positives","n_ag_tests","n_naat_tests","days_lost_waiting")
    return(out)
    
  }
  
  #(C) NAAT Confirmation for Symptomatic Antigen-Negative (Sx/Ag-neg) and Asymptomatic Antigen-Positive (Asx/Ag-pos) Results - each person receives an antigen test and NAAT is used to confirm diagnoses in persons for whom antigen results do not match binary symptom status (e.g., a symptomatic person whose antigen result is negative). 
  naat_sx.agneg_asx.agpos.fun<-function(
    sx.cases = NA,
    sx.noncases  = NA,
    ag.sens.sx = NA,
    ag.sens.asx = NA,
    ag.spec.sx = NA,
    ag.spec.asx = NA,
    naat.sens = NA,
    naat.spec  = NA,
    ag.sens.repeat = NA,
    ag.spec.repeat = NA,
    prev = NA,
    asx.noncases.contacts = NA,
    naat.delay = NA
    
  ){
    detected_cases = 100000*prev*( sx.cases*ag.sens.sx + sx.cases*(1-ag.sens.sx)*naat.sens + (1-sx.cases)*ag.sens.asx*naat.sens )
    missed_cases = 100000*prev*( sx.cases*(1-ag.sens.sx)*(1-naat.sens) + (1-sx.cases)*ag.sens.asx*(1-naat.sens) + (1-sx.cases)*(1-ag.sens.asx) )
    false_positives= 100000*(1-prev)*( sx.noncases*(1-ag.spec.sx) + sx.noncases*ag.spec.sx*(1-naat.spec ) + (1-sx.noncases)*(1-ag.spec.asx)*(1-naat.spec) )
    n_ag_tests = 100000
    n_naat_tests = 100000*( prev*(sx.cases*(1-ag.sens.sx) + (1-sx.cases)*ag.sens.asx) + (1-prev)*(sx.noncases*ag.spec.sx + (1-sx.noncases)*(1-ag.spec.asx) ) )
    days_lost_waiting = 100000*(1-prev)*(sx.noncases*ag.spec.sx +  (1-sx.noncases)*asx.noncases.contacts*(1-ag.spec.asx) + (1-sx.noncases)*(1-asx.noncases.contacts)*(1-ag.spec.asx) )*naat.delay
    
    out=c(detected_cases,missed_cases,false_positives,n_ag_tests,n_naat_tests,days_lost_waiting)
    names(out)<-c("detected_cases","missed_cases","false_positives","n_ag_tests","n_naat_tests","days_lost_waiting")
    return(out)
    
  }
  
  #(D) NAAT Confirmation of Negative Antigen Results (Ag-neg) - each person receives an antigen test and NAAT is used to confirm negative antigen test results.
  naat_agneg.fun<-function(
    sx.cases = NA,
    sx.noncases  = NA,
    ag.sens.sx = NA,
    ag.sens.asx = NA,
    ag.spec.sx = NA,
    ag.spec.asx = NA,
    naat.sens = NA,
    naat.spec  = NA,
    ag.sens.repeat = NA,
    ag.spec.repeat = NA,
    prev = NA,
    asx.noncases.contacts = NA,
    naat.delay = NA
    
  ){
    detected_cases = 100000*prev*( sx.cases*ag.sens.sx + sx.cases*(1-ag.sens.sx)*naat.sens + (1-sx.cases)*ag.sens.asx + (1-sx.cases)*(1-ag.sens.asx)*naat.sens )
    missed_cases = 100000*prev*( sx.cases*(1-ag.sens.sx)*(1-naat.sens) + (1-sx.cases)*(1-ag.sens.asx)*(1-naat.sens) )
    false_positives= 100000*(1- prev)*( sx.noncases*(1-ag.spec.sx) + sx.noncases*ag.spec.sx*(1-naat.spec) + (1-sx.noncases)*(1-ag.spec.asx) + (1-sx.noncases)*ag.spec.asx*(1-naat.spec) )
    n_ag_tests = 100000
    n_naat_tests = 100000*( prev*(sx.cases*(1-ag.sens.sx) +  (1-sx.cases)*(1-ag.sens.asx)) + (1-prev)*(sx.noncases*ag.spec.sx + (1-sx.noncases)*ag.spec.asx ) )
    days_lost_waiting = 100000*(1-prev)*( sx.noncases*ag.spec.sx + (1-sx.noncases)*asx.noncases.contacts*ag.spec.asx )*naat.delay
    
    out=c(detected_cases,missed_cases,false_positives,n_ag_tests,n_naat_tests,days_lost_waiting)
    names(out)<-c("detected_cases","missed_cases","false_positives","n_ag_tests","n_naat_tests","days_lost_waiting")
    return(out)
    
  }
  
  #(E) Repeat Antigen Confirmation of (Ag-neg) - each person receives an antigen test and, for those with initial negative results, a repeat antigen test (performed within approximately 30 minutes of the initial test) is used to confirm negative diagnoses
  repeat_ag_agneg.fun<-function(
    sx.cases = NA,
    sx.noncases  = NA,
    ag.sens.sx = NA,
    ag.sens.asx = NA,
    ag.spec.sx = NA,
    ag.spec.asx = NA,
    naat.sens = NA,
    naat.spec  = NA,
    ag.sens.repeat = NA,
    ag.spec.repeat = NA,
    prev = NA,
    asx.noncases.contacts = NA,
    naat.delay = NA
    
  ){
    detected_cases = 100000*prev*( sx.cases*ag.sens.sx + sx.cases*(1-ag.sens.sx)*ag.sens.repeat + (1-sx.cases)*ag.sens.asx + (1-sx.cases)*(1-ag.sens.asx)*ag.sens.repeat )
    missed_cases = 100000*prev*( sx.cases*(1-ag.sens.sx)*(1- ag.sens.repeat) + (1-sx.cases)*(1-ag.sens.asx)*(1-ag.sens.repeat) )
    false_positives= 100000*(1-prev)*( sx.noncases*(1-ag.spec.sx) + sx.noncases*ag.spec.sx*(1-ag.spec.repeat) + (1-sx.noncases)*(1-ag.spec.asx) + (1-sx.noncases)*ag.spec.asx*(1-ag.spec.repeat) )
    n_ag_tests = 100000 + 100000*( prev*( sx.cases*(1-ag.sens.sx) +  (1-sx.cases)*(1-ag.sens.asx) )   + (1- prev)*(sx.noncases*ag.spec.sx +(1-sx.noncases)*ag.spec.asx ) )
    n_naat_tests = 0
    days_lost_waiting = 0
    
    out=c(detected_cases,missed_cases,false_positives,n_ag_tests,n_naat_tests,days_lost_waiting)
    names(out)<-c("detected_cases","missed_cases","false_positives","n_ag_tests","n_naat_tests","days_lost_waiting")
    return(out)
    
    
  }
  
  #(F) NAAT for Asymptomatic Persons (Asx) & Symptomatic Persons with Positive Antigen Results (Sx/Ag-pos) - asymptomatic persons receive a NAAT; symptomatic persons receive an antigen test followed by a NAAT for those with positive antigen results
  naat_asx_sx.agpos.fun<-function(
    sx.cases = NA,
    sx.noncases  = NA,
    ag.sens.sx = NA,
    ag.sens.asx = NA,
    ag.spec.sx = NA,
    ag.spec.asx = NA,
    naat.sens = NA,
    naat.spec  = NA,
    ag.sens.repeat = NA,
    ag.spec.repeat = NA,
    prev = NA,
    asx.noncases.contacts = NA,
    naat.delay = NA
    
  ){
    detected_cases = 100000*prev*( sx.cases*ag.sens.sx*naat.sens + (1-sx.cases)*naat.sens )
    missed_cases = 100000*prev*( sx.cases*(1-ag.sens.sx) + sx.cases*ag.sens.sx*(1-naat.sens) + (1-sx.cases)*(1-naat.sens) )
    false_positives= 100000*(1- prev)*( sx.noncases*(1-ag.spec.sx)*(1-naat.spec) + (1-sx.noncases)*(1-naat.spec) )
    n_ag_tests = 100000*( prev*sx.cases  + (1-prev)*sx.noncases  )
    n_naat_tests = 100000*( prev*(sx.cases*ag.sens.sx + (1-sx.cases) ) + (1-prev)*(sx.noncases*(1- ag.spec.sx) + (1-sx.noncases) ) )
    days_lost_waiting = 100000*(1-prev)*(sx.noncases*(1-ag.spec.sx) + (1-sx.noncases)*asx.noncases.contacts )*naat.delay
    
    out=c(detected_cases,missed_cases,false_positives,n_ag_tests,n_naat_tests,days_lost_waiting)
    names(out)<-c("detected_cases","missed_cases","false_positives","n_ag_tests","n_naat_tests","days_lost_waiting")
    return(out)
    
  }
