## ---------------------------------------------
## 1. Initial conditions
## ---------------------------------------------

initial.state.y <- c(
  S.y      = 9990,
  E.y      = 10,
  Ia.sh.y  = 0,
  Im.syU.y = 0,
  Im.sh.y  = 0,
  Im.syT.y = 0,
  Im.abx.y = 0,
  Is.syU.y = 0,
  Is.sh.y  = 0,
  Is.syT.y = 0,
  Is.abx.y = 0,
  R.y      = 0,
  D.y      = 0,
  R.abx.y  = 0
)

initial.state.o <- c(
  S.o      = 9990,
  E.o      = 10,
  Ia.sh.o  = 0,
  Im.syU.o = 0,
  Im.sh.o  = 0,
  Im.syT.o = 0,
  Im.abx.o = 0,
  Is.syU.o = 0,
  Is.sh.o  = 0,
  Is.syT.o = 0,
  Is.abx.o = 0,
  R.o      = 0,
  D.o      = 0,
  R.abx.o  = 0
)

## ---------------------------------------------
## 2. Parameter values
##    (THESE ARE ALL PLACEHODLER PARAMETERS)
## ---------------------------------------------

## shared-ish choices
sigma      <- 1/2       # 2-day latent period
beta.y     <- 1.5      
beta.o     <- 1.2

v.a.y      <- 0.3
v.m.y      <- 1.0
v.sh.y     <- 2.0
v.abx.y    <- 0.5

epsilon.a.y   <- 0.3
epsilon.s.y   <- 0.3
epsilon.m.T.y <- 0.5
epsilon.s.T.y <- 0.7

alpha.m.y  <- 1/3
alpha.s.y  <- 1/3
tau.y      <- 1/4
q.y        <- 0.6
delta.y    <- 0.8

theta.m.y  <- 0.05
theta.s.y  <- 0.2

gamma.a.y      <- 1/5
gamma.m.y      <- 1/7
gamma.s.y      <- 1/10
gamma.m.abx.y  <- 1/4
gamma.s.abx.y  <- 1/6

mu.m.y     <- 1/30
mu.s.y     <- 1/15

## adults – here just making them slightly different
v.a.o      <- 0.25
v.m.o      <- 0.9
v.sh.o     <- 1.8
v.abx.o    <- 0.4

epsilon.a.o   <- 0.25
epsilon.s.o   <- 0.35
epsilon.m.T.o <- 0.6
epsilon.s.T.o <- 0.75

alpha.m.o  <- 1/3
alpha.s.o  <- 1/3
tau.o      <- 1/4
q.o        <- 0.5
delta.o    <- 0.8

theta.m.o  <- 0.06
theta.s.o  <- 0.22

gamma.a.o      <- 1/5
gamma.m.o      <- 1/8
gamma.s.o      <- 1/11
gamma.m.abx.o  <- 1/4
gamma.s.abx.o  <- 1/6

mu.m.o     <- 1/28
mu.s.o     <- 1/14

omega <- 1/365  
l     <- 0.6  

## ---------------------------------------------
## 3. Call the model
##    NOTE: sigma appears twice in your function signature.
##    Here we pass it once by name (y) and once positionally (o).
## ---------------------------------------------

set.seed(123)

source("AgeStrAll.R")

age.out <- AgeStructSEIRfunct(
  ## youth block
  beta.y          = beta.y,
  sigma           = sigma,
  v.a.y           = v.a.y,
  v.m.y           = v.m.y,
  v.sh.y          = v.sh.y,
  v.abx.y         = v.abx.y,
  epsilon.a.y     = epsilon.a.y,
  epsilon.s.y     = epsilon.s.y,
  epsilon.m.T.y   = epsilon.m.T.y,
  epsilon.s.T.y   = epsilon.s.T.y,
  alpha.m.y       = alpha.m.y,
  alpha.s.y       = alpha.s.y,
  tau.y           = tau.y,
  q.y             = q.y,
  delta.y         = delta.y,
  theta.m.y       = theta.m.y,
  theta.s.y       = theta.s.y,
  gamma.a.y       = gamma.a.y,
  gamma.m.y       = gamma.m.y,
  gamma.s.y       = gamma.s.y,
  gamma.m.abx.y   = gamma.m.abx.y,
  gamma.s.abx.y   = gamma.s.abx.y,
  mu.m.y          = mu.m.y,
  mu.s.y          = mu.s.y,
  initial.state.y = initial.state.y,
  
  ## adult block – after initial.state.y the args are positional
  beta.o          = beta.o,
  v.a.o           = v.a.o,
  v.m.o           = v.m.o,
  v.sh.o          = v.sh.o,
  v.abx.o         = v.abx.o,
  epsilon.a.o     = epsilon.a.o,
  epsilon.s.o     = epsilon.s.o,
  epsilon.m.T.o   = epsilon.m.T.o,
  epsilon.s.T.o   = epsilon.s.T.o,
  alpha.m.o       = alpha.m.o,
  alpha.s.o       = alpha.s.o,
  tau.o           = tau.o,
  q.o             = q.o,
  delta.o         = delta.o,
  theta.m.o       = theta.m.o,
  theta.s.o       = theta.s.o,
  gamma.a.o       = gamma.a.o,
  gamma.m.o       = gamma.m.o,
  gamma.s.o       = gamma.s.o,
  gamma.m.abx.o   = gamma.m.abx.o,
  gamma.s.abx.o   = gamma.s.abx.o,
  mu.m.o          = mu.m.o,
  mu.s.o          = mu.s.o,
  initial.state.o = initial.state.o,
  omega           = omega,
  l               = l,
  
  ## options
  step.size       = 1,
  freq.dependent  = TRUE
)

## ---------------------------------------------
## 4. Inspect output
## ---------------------------------------------

str(age.out)
head(age.out$SEIR.output)
head(age.out$incidences)

df <- age.out$SEIR.output %>%
  pivot_longer(-time, names_to = "compartment", values_to = "value")

ggplot(df, aes(x = time, y = value)) +
  geom_line() +
  facet_wrap(~ compartment, scales = "free_y") +
  theme_bw() +
  labs(title = "SEIR Compartments Over Time",
       x = "Time",
       y = "Count")
