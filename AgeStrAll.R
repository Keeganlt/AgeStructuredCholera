library(dplyr)
library(ggplot2)
library(lhs)
library(dplyr) #need for bind_rows for incidences_sum
library(reshape2) #need for melt for total.inc_long
library(tidyr)

# runSEIR takes the SEIR model parameters and initial condition 
# and returns a time series for each group
AgeStructSEIRfunct <- function(beta.y, sigma, v.a.y, v.m.y, v.sh.y, v.abx.y,
                       epsilon.a.y, epsilon.s.y, epsilon.m.T.y, epsilon.s.T.y,
                       alpha.m.y, alpha.s.y,
                       tau.y, q.y, delta.y,
                       theta.m.y, theta.s.y,
                       gamma.a.y, gamma.m.y, gamma.s.y,
                       gamma.m.abx.y, gamma.s.abx.y,
                       mu.m.y, mu.s.y,
                       initial.state.y,
                       beta.o, v.a.o, v.m.o, v.sh.o, v.abx.o,
                       epsilon.a.o, epsilon.s.o, epsilon.m.T.o, epsilon.s.T.o,
                       alpha.m.o, alpha.s.o,
                       tau.o, q.o, delta.o,
                       theta.m.o, theta.s.o,
                       gamma.a.o, gamma.m.o, gamma.s.o,
                       gamma.m.abx.o, gamma.s.abx.o,
                       mu.m.o, mu.s.o,
                       initial.state.o,
                       omega,
                       l,
                       step.size = 1,
                       freq.dependent=TRUE){

  #If the model is frequency dependent we modify beta based on the total population size
  beta.divisor.y <- ifelse(freq.dependent,
                           initial.state.y["S.y"]+
                             initial.state.y["E.y"]+
                             initial.state.y["Ia.sh.y"]+
                             initial.state.y["Im.syU.y"]+initial.state.y["Im.sh.y"]+
                             initial.state.y["Im.syT.y"]+initial.state.y["Im.abx.y"]+
                             initial.state.y["Is.syU.y"]+initial.state.y["Is.sh.y"]+
                             initial.state.y["Is.syT.y"]+initial.state.y["Is.abx.y"]+
                             initial.state.y["R.y"]+initial.state.y["R.abx.y"]+initial.state.y["D.y"],
                           1)
  #If the model is frequency dependent we modify beta based on the total population size
  beta.divisor.o <- ifelse(freq.dependent,
                           initial.state.o["S.o"]+
                             initial.state.o["E.o"]+
                             initial.state.o["Ia.sh.o"]+
                             initial.state.o["Im.syU.o"]+initial.state.o["Im.sh.o"]+
                             initial.state.o["Im.syT.o"]+initial.state.o["Im.abx.o"]+
                             initial.state.o["Is.syU.o"]+initial.state.o["Is.sh.o"]+
                             initial.state.o["Is.syT.o"]+initial.state.o["Is.abx.o"]+
                             initial.state.o["R.o"]+initial.state.o["R.abx.o"]+initial.state.o["D.o"],
                           1)
  
  #create the parameter vector.
  param <- c(beta.scaled.y = beta.y/beta.divisor.y, 
             sigma = sigma, 
             v.a.y = v.a.y, 
             v.m.y = v.m.y, 
             v.sh.y = v.sh.y, 
             v.abx.y = v.abx.y,
             epsilon.a.y = epsilon.a.y, 
             epsilon.s.y = epsilon.s.y, 
             epsilon.m.T.y = epsilon.m.T.y, 
             epsilon.s.T.y = epsilon.s.T.y,
             alpha.m.y = alpha.m.y, 
             alpha.s.y = alpha.s.y,
             tau.y = tau.y,
             q.y = q.y,
             delta.y = delta.y,
             theta.m.y = theta.m.y, 
             theta.s.y = theta.s.y,
             gamma.a.y = gamma.a.y, 
             gamma.m.y = gamma.m.y, 
             gamma.s.y = gamma.s.y,
             gamma.m.abx.y = gamma.m.abx.y, 
             gamma.s.abx.y = gamma.s.abx.y,
             mu.m.y = mu.m.y, 
             mu.s.y = mu.s.y,
             beta.scaled.o = beta.o/beta.divisor.o, 
             sigma = sigma, 
             v.a.o = v.a.o, 
             v.m.o = v.m.o, 
             v.sh.o = v.sh.o, 
             v.abx.o = v.abx.o,
             epsilon.a.o = epsilon.a.o, 
             epsilon.s.o = epsilon.s.o, 
             epsilon.m.T.o = epsilon.m.T.o, 
             epsilon.s.T.o = epsilon.s.T.o,
             alpha.m.o = alpha.m.o, 
             alpha.s.o = alpha.s.o,
             tau.o = tau.o,
             q.o = q.o,
             delta.o = delta.o,
             theta.m.o = theta.m.o, 
             theta.s.o = theta.s.o,
             gamma.a.o = gamma.a.o, 
             gamma.m.o = gamma.m.o, 
             gamma.s.o = gamma.s.o,
             gamma.m.abx.o = gamma.m.abx.o, 
             gamma.s.abx.o = gamma.s.abx.o,
             mu.m.o = mu.m.o, 
             mu.s.o = mu.s.o,
             omega = omega,  
             l = l)
  
  #Since we are not using a fancy solver we will need to run this on our own.
  #note that the simulation ends once there are 0 people in groups E or I
  t <- 0

  y.y <- initial.state.y

  y.o <- initial.state.o
  SEIR.output <- matrix(ncol=length(initial.state.o) + length(initial.state.y)+1, nrow=1)
  #colnames(SEIR.output)
  name.s <- c("time",
                             "S.y", 
                             "E.y",
                             "Ia.sh.y",
                             "Im.syU.y", 
                             "Im.sh.y", 
                             "Im.syT.y", 
                             "Im.abx.y",
                             "Is.syU.y", 
                             "Is.sh.y", 
                             "Is.syT.y", 
                             "Is.abx.y",
                             "R.y", 
                             "D.y", 
                             "R.abx.y", 
                             "S.o", 
                             "E.o",
                             "Ia.sh.o",
                             "Im.syU.o", 
                             "Im.sh.o", 
                             "Im.syT.o", 
                             "Im.abx.o",
                             "Is.syU.o", 
                             "Is.sh.o", 
                             "Is.syT.o", 
                             "Is.abx.o",
                             "R.o", 
                             "D.o", 
                             "R.abx.o")
  colnames(SEIR.output) <-name.s
  SEIR.output[1,] <- c(t, y.y, y.o)
  incidences <- matrix(ncol=18+18+1, nrow=1)
  colnames(incidences) <- c("time",
                            "incident.exposed.y",
                            "incident.cases.a.y",
                            "incident.cases.symp.mU.y", 
                            "incident.cases.symp.mT.y",
                            "incident.cases.symp.sU.y", 
                            "incident.cases.symp.sT.y",
                            "shedding.cases.mU.y", 
                            "shedding.cases.mT.y", 
                            "abx.cases.m.y",
                            "shedding.cases.s.y", 
                            "abx.cases.s.y",
                            "nat.recovered.a.y",
                            "nat.recovered.m.y", 
                            "died.m.y", 
                            "abx.recovered.m.y",
                            "nat.recovered.s.y", 
                            "died.s.y", 
                            "abx.recovered.s.y",
                            "incident.exposed.o",
                            "incident.cases.a.o",
                            "incident.cases.symp.mU.o", 
                            "incident.cases.symp.mT.o",
                            "incident.cases.symp.sU.o", 
                            "incident.cases.symp.sT.o",
                            "shedding.cases.mU.o", 
                            "shedding.cases.mT.o", 
                            "abx.cases.m.o",
                            "shedding.cases.s.o", 
                            "abx.cases.s.o",
                            "nat.recovered.a.o",
                            "nat.recovered.m.o", 
                            "died.m.o", 
                            "abx.recovered.m.o",
                            "nat.recovered.s.o", 
                            "died.s.o", 
                            "abx.recovered.s.o")
  incidences[1,] <- rep(0,18+18+1)
  
  while(
        y.y["E.y"]>0 || 
    y.y["Ia.sh.y"]>0 ||
    y.y["Im.syU.y"]>0 || 
    y.y["Im.sh.y"]>0 || 
    y.y["Im.syT.y"]>0 || 
    y.y["Im.abx.y"]>0 ||
    y.y["Is.syU.y"]>0 || 
    y.y["Is.sh.y"]>0 || 
    y.y["Is.syT.y"]>0 || 
    y.y["Is.abx.y"]>0 ||
    y.o["E.o"]>0 || 
    y.o["Ia.sh.o"]>0 ||
    y.o["Im.syU.o"]>0 || 
    y.o["Im.sh.o"]>0 || 
    y.o["Im.syT.o"]>0 || 
    y.o["Im.abx.o"]>0 ||
    y.o["Is.syU.o"]>0 || 
    y.o["Is.sh.o"]>0 || 
    y.o["Is.syT.o"]>0 || 
    y.o["Is.abx.o"]>0) {
    
    t <- t+step.size 
    
    lambda <- (param["beta.scaled.y"]*(y.y["Is.syU.y"] + y.y["Is.syT.y"])) + 
      (v.sh.y*param["beta.scaled.y"]*y.y["Is.sh.y"]) + 
      (v.a.y*param["beta.scaled.y"]*y.y["Ia.sh.y"]) + 
      ((v.m.y*param["beta.scaled.y"])*(y.y["Im.syU.y"]+y.y["Im.syT.y"])) + 
      (v.sh.y*v.m.y*param["beta.scaled.y"]*y.y["Im.sh.y"]) + 
      (v.m.y*v.abx.y*param["beta.scaled.y"]*y.y["Im.abx.y"]) +
      (param["beta.scaled.o"]*(y.o["Is.syU.o"] + y.o["Is.syT.o"])) + 
      (v.sh.o*param["beta.scaled.o"]*y.o["Is.sh.o"]) + 
      (v.a.o*param["beta.scaled.o"]*y.o["Ia.sh.o"]) + 
      ((v.m.o*param["beta.scaled.o"])*(y.o["Im.syU.o"]+y.o["Im.syT.o"])) + 
      (v.sh.o*v.m.o*param["beta.scaled.o"]*y.o["Im.sh.o"]) + 
      (v.m.o*v.abx.o*param["beta.scaled.o"]*y.o["Im.abx.o"])
    
    lambda
    
    #at this point, lambda is basically beta, looks complicated since is so many different infected types
    #is a number, should look like a number. print to look at
    
    #calculate the probability of infection and recovery in this time step
    #probability that thing happens in a given time step, updates every time step. range 0-1
    #need one of these per arrow, moderate/severe specific etc
## -------- youth 
    Pr.age.y <- 1-exp(-step.size * param["omega"])
    
    Pr.expose.y <- 1-exp(-step.size*lambda)
    
    Pr.infect.all.y <- 1-exp(-step.size*param["sigma"])
    
    Pr.mU.shed.y <- 1-exp(-step.size*param["alpha.m.y"]) #the rate at which moderates who don't seek care move to shed
    Pr.mT.shed.y <- 1-exp(-step.size*param["alpha.m.y"]) #this is only rate of move to shed
    Pr.mT.abx.y <- 1-exp(-step.size*param["delta.y"]*param["tau.y"]) #this is only rate of moderates who seek care move to abx
    
    Pr.sU.shed.y <- 1-exp(-step.size*param["alpha.s.y"]) #the rate at which severes who don't seek care move to shed
    Pr.sT.abx.y <- 1-exp(-step.size*param["tau.y"]) #the rate at which severes who do seek care move to abx
    
    Pr.a.recover.y <- 1-exp(-step.size*param["gamma.a.y"]) #rate at which asymptomatics recover
    Pr.m.shed.recover.y <- 1-exp(-step.size*param["gamma.m.y"]) #this only incorporates the rate at which shedding moderates recover
    Pr.m.shed.die.y <- 1-exp(-step.size*param["mu.m.y"]) #this only incorporates the rate at which shedding moderates die
    
    Pr.m.abx.recover.y <- 1-exp(-step.size*param["gamma.m.abx.y"]) #rate at which moderates who receive abx recover
    Pr.s.shed.recover.y <- 1-exp(-step.size*param["gamma.s.y"]) #this only incorporates the rate of recover among shedding severes
    Pr.s.shed.die.y <- 1-exp(-step.size*param["mu.s.y"]) #this only incorporates the rate of die among shedding severes
    
    Pr.s.abx.recover.y <- 1-exp(-step.size*param["gamma.s.abx.y"]) #rate of recover among severes who receive abx
    
    #draw random variable from binomial distribution for new number of
    #using prob from above, is the actual number based on prob and denom (num ppl in that box)
    #should be whole numbers. if aren't something's wrong
    #rbinom(num observations, num trials, prob success each trial)
    incident.exposed.y <- rbinom(1, y.y["S.y"], Pr.expose.y) #the number who go S to E
    S.y.remaining <- y.y["S.y"] - incident.exposed.y
    aged.y <- rbinom(1, S.y.remaining, Pr.age.y)
    
    incident.cases.all.y <- rbinom(1, y.y["E.y"], Pr.infect.all.y) #this is the number that leave E
    incident.cases.a.y <- round(incident.cases.all.y*param["epsilon.a.y"],0) #this the number that enter I_a
    
    #now only have (number that leave)-(number that go to A) left to go to symptomatics
    incident.cases.symp.y <- (incident.cases.all.y-incident.cases.a.y) #this is the number that go to symptomatic I's
    incident.cases.symp.s.y <- round(incident.cases.symp.y*param["epsilon.s.y"],0) #this is the number that go to severe symptomatic I's
    incident.cases.symp.m.y <- (incident.cases.symp.y-incident.cases.symp.s.y) #this is the number that goes to moderate symptomatic I's
    incident.cases.symp.mT.y <- round(incident.cases.symp.m.y*param["epsilon.m.T.y"],0) #this is the number that go to moderate seek care
    incident.cases.symp.mU.y <- (incident.cases.symp.m.y-incident.cases.symp.mT.y) #this is the number that go to moderate don't seek care
    incident.cases.symp.sT.y <- round(incident.cases.symp.s.y*param["epsilon.s.T.y"],0) #this is the number that go to severe seek care
    incident.cases.symp.sU.y <- (incident.cases.symp.s.y-incident.cases.symp.sT.y) #this is the number that go to severe don't seek care
    

    
    shedding.cases.mU.y <- rbinom(1, y.y["Im.syU.y"], Pr.mU.shed.y) #the number who move to moderate shedding
    m.who.will.abx.y <- round(y.y["Im.syT.y"]*param["q.y"],0) #the number of moderates who seek care who will eventually receive abx#are treating q like a proportion, though this not actually correct. the math works though and conceptually is ok via prop.m.abx equation
    shedding.cases.mT.y <- rbinom(1, y.y["Im.syT.y"]-m.who.will.abx.y, Pr.mT.shed.y) #the number of moderates who seek care who don't receive abx and move to shed
    abx.cases.m.y <- rbinom(1, m.who.will.abx.y, Pr.mT.abx.y) #the number of moderates who receive abx at this step
    
    
    shedding.cases.s.y <- rbinom(1, y.y["Is.syU.y"], Pr.sU.shed.y) #the number of severes who don't seek care and move to shed
    abx.cases.s.y <- rbinom(1, y.y["Is.syT.y"], Pr.sT.abx.y) #the number of severes who seek care and receive abx
    
    nat.recovered.a.y <- rbinom(1, y.y["Ia.sh.y"], Pr.a.recover.y) #the number of asymptomatics who recover
    m.who.will.die.y <- round(y.y["Im.sh.y"]*param["theta.m.y"],0) #this is the whole number of moderates who will eventually die
    nat.recovered.m.y <- rbinom(1, y.y["Im.sh.y"]-m.who.will.die.y, Pr.m.shed.recover.y) #this is the actual number of I_sh who naturally recover. this incorporates the stochastic prob that will actually move in this step #rbinom requires whole numbers
    died.m.y <- rbinom(1, m.who.will.die.y, Pr.m.shed.die.y) #the number of moderates who move from I_sh to dead
    
    abx.recovered.m.y <- rbinom(1, y.y["Im.abx.y"], Pr.m.abx.recover.y) #number of moderates who move from abx to recover in this step
    s.who.will.die.y <- round(y.y["Is.sh.y"]*param["theta.s.y"],0) #this is the whole number of severes who will die
    nat.recovered.s.y <- rbinom(1, y.y["Is.sh.y"]-s.who.will.die.y, Pr.s.shed.recover.y) #number of severes who move from Is_sh to recover w/o abx
    died.s.y <- rbinom(1, s.who.will.die.y, Pr.s.shed.die.y) #number of Is_sh who die in this time step
    abx.recovered.s.y <- rbinom(1, y.y["Is.abx.y"], Pr.s.abx.recover.y) #number of severes who recover after abx in this time step


## -------- Adults


    Pr.expose.o <- 1-exp(-step.size * param["l"] * lambda)
    
    Pr.infect.all.o <- 1-exp(-step.size*param["sigma"])
    
    Pr.mU.shed.o <- 1-exp(-step.size*param["alpha.m.o"]) #the rate at which moderates who don't seek care move to shed
    Pr.mT.shed.o <- 1-exp(-step.size*param["alpha.m.o"]) #this is only rate of move to shed
    Pr.mT.abx.o <- 1-exp(-step.size*param["delta.o"]*param["tau.o"]) #this is only rate of moderates who seek care move to abx
    
    Pr.sU.shed.o <- 1-exp(-step.size*param["alpha.s.o"]) #the rate at which severes who don't seek care move to shed
    Pr.sT.abx.o <- 1-exp(-step.size*param["tau.o"]) #the rate at which severes who do seek care move to abx
    
    Pr.a.recover.o <- 1-exp(-step.size*param["gamma.a.o"]) #rate at which asymptomatics recover
    Pr.m.shed.recover.o <- 1-exp(-step.size*param["gamma.m.o"]) #this only incorporates the rate at which shedding moderates recover
    Pr.m.shed.die.o <- 1-exp(-step.size*param["mu.m.o"]) #this only incorporates the rate at which shedding moderates die
    
    Pr.m.abx.recover.o <- 1-exp(-step.size*param["gamma.m.abx.o"]) #rate at which moderates who receive abx recover
    Pr.s.shed.recover.o <- 1-exp(-step.size*param["gamma.s.o"]) #this only incorporates the rate of recover among shedding severes
    Pr.s.shed.die.o <- 1-exp(-step.size*param["mu.s.o"]) #this only incorporates the rate of die among shedding severes
    
    Pr.s.abx.recover.o <- 1-exp(-step.size*param["gamma.s.abx.o"]) #rate of recover among severes who receive abx
    
    #draw random variable from binomial distribution for new number of
    #using prob from above, is the actual number based on prob and denom (num ppl in that box)
    #should be whole numbers. if aren't something's wrong
    #rbinom(num observations, num trials, prob success each trial)
    incident.exposed.o <- rbinom(1, y.o["S.o"], Pr.expose.o) #the number who go S to E
    
    incident.cases.all.o <- rbinom(1, y.o["E.o"], Pr.infect.all.o) #this is the number that leave E
    incident.cases.a.o <- round(incident.cases.all.o*param["epsilon.a.o"],0) #this the number that enter I_a
    
    #now only have (number that leave)-(number that go to A) left to go to symptomatics
    incident.cases.symp.o <- (incident.cases.all.o-incident.cases.a.o) #this is the number that go to symptomatic I's
    incident.cases.symp.s.o <- round(incident.cases.symp.o*param["epsilon.s.o"],0) #this is the number that go to severe symptomatic I's
    incident.cases.symp.m.o <- (incident.cases.symp.o-incident.cases.symp.s.o) #this is the number that goes to moderate symptomatic I's
    incident.cases.symp.mT.o <- round(incident.cases.symp.m.o*param["epsilon.m.T.o"],0) #this is the number that go to moderate seek care
    incident.cases.symp.mU.o <- (incident.cases.symp.m.o-incident.cases.symp.mT.o) #this is the number that go to moderate don't seek care
    incident.cases.symp.sT.o <- round(incident.cases.symp.s.o*param["epsilon.s.T.o"],0) #this is the number that go to severe seek care
    incident.cases.symp.sU.o <- (incident.cases.symp.s.o-incident.cases.symp.sT.o) #this is the number that go to severe don't seek care
    

    
    shedding.cases.mU.o <- rbinom(1, y.o["Im.syU.o"], Pr.mU.shed.o) #the number who move to moderate shedding
    m.who.will.abx.o <- round(y.o["Im.syT.o"]*param["q.o"],0) #the number of moderates who seek care who will eventually receive abx#are treating q like a proportion, though this not actually correct. the math works though and conceptually is ok via prop.m.abx equation
    shedding.cases.mT.o <- rbinom(1, y.o["Im.syT.o"]-m.who.will.abx.o, Pr.mT.shed.o) #the number of moderates who seek care who don't receive abx and move to shed
    abx.cases.m.o <- rbinom(1, m.who.will.abx.o, Pr.mT.abx.o) #the number of moderates who receive abx at this step
    
    
    shedding.cases.s.o <- rbinom(1, y.o["Is.syU.o"], Pr.sU.shed.o) #the number of severes who don't seek care and move to shed
    abx.cases.s.o <- rbinom(1, y.o["Is.syT.o"], Pr.sT.abx.o) #the number of severes who seek care and receive abx
    
    nat.recovered.a.o <- rbinom(1, y.o["Ia.sh.o"], Pr.a.recover.o) #the number of asymptomatics who recover
    m.who.will.die.o <- round(y.o["Im.sh.o"]*param["theta.m.o"],0) #this is the whole number of moderates who will eventually die
    nat.recovered.m.o <- rbinom(1, y.o["Im.sh.o"]-m.who.will.die.o, Pr.m.shed.recover.o) #this is the actual number of I_sh who naturally recover. this incorporates the stochastic prob that will actually move in this step #rbinom requires whole numbers
    died.m.o <- rbinom(1, m.who.will.die.o, Pr.m.shed.die.o) #the number of moderates who move from I_sh to dead
    
    abx.recovered.m.o <- rbinom(1, y.o["Im.abx.o"], Pr.m.abx.recover.o) #number of moderates who move from abx to recover in this step
    s.who.will.die.o <- round(y.o["Is.sh.o"]*param["theta.s.o"],0) #this is the whole number of severes who will die
    nat.recovered.s.o <- rbinom(1, y.o["Is.sh.o"]-s.who.will.die.o, Pr.s.shed.recover.o) #number of severes who move from Is_sh to recover w/o abx
    died.s.o <- rbinom(1, s.who.will.die.o, Pr.s.shed.die.o) #number of Is_sh who die in this time step
    abx.recovered.s.o <- rbinom(1, y.o["Is.abx.o"], Pr.s.abx.recover.o) #number of severes who recover after abx in this time step
    
    #Find the deltas for each compartment

    dS.y <- -incident.exposed.y - aged.y
    dE.y <- incident.exposed.y - incident.cases.a.y - 
      incident.cases.symp.mU.y - incident.cases.symp.mT.y - 
      incident.cases.symp.sU.y - incident.cases.symp.sT.y
    dIa.sh.y <- incident.cases.a.y - nat.recovered.a.y
    

    dIm.syU.y <- incident.cases.symp.mU.y - shedding.cases.mU.y
    dIm.sh.y <- shedding.cases.mU.y + shedding.cases.mT.y - nat.recovered.m.y - died.m.y
    dIm.syT.y <- incident.cases.symp.mT.y - shedding.cases.mT.y - abx.cases.m.y
    dIm.abx.y <- abx.cases.m.y - abx.recovered.m.y
    
    dIs.syU.y <- incident.cases.symp.sU.y - shedding.cases.s.y
    dIs.sh.y <- shedding.cases.s.y - nat.recovered.s.y - died.s.y
    dIs.syT.y <- incident.cases.symp.sT.y - abx.cases.s.y
    dIs.abx.y <- abx.cases.s.y - abx.recovered.s.y
    
    dR.y <- nat.recovered.a.y + nat.recovered.m.y + nat.recovered.s.y
    dD.y <- died.m.y + died.s.y
    dR.abx.y <- abx.recovered.m.y + abx.recovered.s.y
    
    sum(dIa.sh.y, dIm.syU.y, dIm.syT.y, dIs.syU.y, dIs.syT.y)


    dS.o <- -incident.exposed.o + aged.y
    dE.o <- incident.exposed.o - incident.cases.a.o - 
      incident.cases.symp.mU.o - incident.cases.symp.mT.o - 
      incident.cases.symp.sU.o - incident.cases.symp.sT.o
    dIa.sh.o <- incident.cases.a.o - nat.recovered.a.o

    
    dIm.syU.o <- incident.cases.symp.mU.o - shedding.cases.mU.o
    dIm.sh.o <- shedding.cases.mU.o + shedding.cases.mT.o - nat.recovered.m.o - died.m.o
    dIm.syT.o <- incident.cases.symp.mT.o - shedding.cases.mT.o - abx.cases.m.o
    dIm.abx.o <- abx.cases.m.o - abx.recovered.m.o
    
    dIs.syU.o <- incident.cases.symp.sU.o - shedding.cases.s.o
    dIs.sh.o <- shedding.cases.s.o - nat.recovered.s.o - died.s.o
    dIs.syT.o <- incident.cases.symp.sT.o - abx.cases.s.o
    dIs.abx.o <- abx.cases.s.o - abx.recovered.s.o
    
    dR.o <- nat.recovered.a.o + nat.recovered.m.o + nat.recovered.s.o
    dD.o <- died.m.o + died.s.o
    dR.abx.o <- abx.recovered.m.o + abx.recovered.s.o
    
    sum(dIa.sh.o, dIm.syU.o, dIm.syT.o, dIs.syU.o, dIs.syT.o)

    deltas.y <- c(dS.y, dE.y, dIa.sh.y, dIm.syU.y, dIm.sh.y, dIm.syT.y, dIm.abx.y,
                  dIs.syU.y, dIs.sh.y, dIs.syT.y, dIs.abx.y, dR.y, dD.y, dR.abx.y) # calculate step sizes. this is net change for each box
    
    
    deltas.o <- c(dS.o, dE.o, dIa.sh.o, dIm.syU.o, dIm.sh.o, dIm.syT.o, dIm.abx.o,
                  dIs.syU.o, dIs.sh.o, dIs.syT.o, dIs.abx.o, dR.o, dD.o, dR.abx.o) # calculate step sizes. this is net change for each box
    
    z.y <- c(incident.exposed.y,
             incident.cases.a.y,
             incident.cases.symp.mU.y, incident.cases.symp.mT.y,
             incident.cases.symp.sU.y, incident.cases.symp.sT.y,
             shedding.cases.mU.y, shedding.cases.mT.y, abx.cases.m.y,
             shedding.cases.s.y, abx.cases.s.y,
             nat.recovered.a.y,
             nat.recovered.m.y, died.m.y, abx.recovered.m.y,
             nat.recovered.s.y, died.s.y, abx.recovered.s.y)
    #this is the incident (new) number of ppl along each arrow at each step
    
    z.o <- c(incident.exposed.o,
             incident.cases.a.o,
             incident.cases.symp.mU.o, incident.cases.symp.mT.o,
             incident.cases.symp.sU.o, incident.cases.symp.sT.o,
             shedding.cases.mU.o, shedding.cases.mT.o, abx.cases.m.o,
             shedding.cases.s.o, abx.cases.s.o,
             nat.recovered.a.o,
             nat.recovered.m.o, died.m.o, abx.recovered.m.o,
             nat.recovered.s.o, died.s.o, abx.recovered.s.o)
    #this is the incident (new) number of ppl along each arrow at each step
    # update states
    y.y <- y.y + deltas.y 
    y.o <- y.o + deltas.o 
    
    # append current state to SEIR.output
    SEIR.output <- rbind(
      SEIR.output,
      c(time = t, y.y, y.o)
    )
    
    # append current incidences
    incidences <- rbind(
      incidences,
      c(time = t, z.y, z.o)
    )
    
  }  # end while
  
  combo <- list(
    SEIR.output = as.data.frame(SEIR.output),
    incidences  = as.data.frame(incidences)
  )
  return(combo)
}
