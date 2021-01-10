library(tidyverse) 

meas <- read_csv(paste0("https://kingaa.github.io/sbied/stochsim/",
                                  "Measles_Consett_1948.csv")) %>% select(week,reports=cases)

meas %>% as.data.frame() %>% head()

plot(meas$week, meas$reports, )

meas %>%
  ggplot(aes(x=week, y=reports))+
  geom_line()+
  geom_point()

#### SIR MODEL ###
sir_step <- function (S, I, R, N, Beta, mu_IR, delta.t, ...) { 
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t)) 
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-mu_IR*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR 
  R <- R + dN_IR
  c(S = S, I = I, R = R)
}

#### Init. Conditions I(0)=1, S(0)= N*eta, R(0)=N*(1-eta), eta is initial fraction of infected individuals
#    as is gonna be estimated ###
sir_rinit <- function (N, eta, ...) 
{
  c(S = round(N*eta), I = 1, R = round(N*(1-eta)))
}


library(pomp)
measSIR <- pomp(data = meas, times="week",t0=0, rprocess=euler(sir_step,delta.t=1/7), rinit=sir_rinit) 

##### MODEL ASSUMING REPORT IS INCIDENCE ##### 
sir_step <- function (S, I, R, H, N, Beta, mu_IR, delta.t, ...) 
{
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t)) 
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-mu_IR*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  H <- H + dN_IR;
  c(S = S, I = I, R = R, H = H)
}

sir_rinit <- function (N, eta, ...) {
  c(S = round(N*eta), I = 1, R = round(N*(1-eta)), H = 0)
}
measSIR <- pomp(data=meas, times="week", t0=0,
  rprocess=euler(sir_step,delta.t=1/7),
  rinit=sir_rinit,accumvars="H" )


#### Report Model is modeled as binomial ####
# Observation measure is Binomial  BIN( H, ) 
sir_dmeas <- function (reports, H, rho, log, ...) { dbinom(x=reports, size=H, prob=rho, log=log) }
sir_rmeas <- function (H, rho, ...) { c(reports=rbinom(n=1, size=H, prob=rho)) }


measSIR = pomp(data=meas, times="week", t0=0,
               rmeasure=sir_rmeas, dmeasure=sir_dmeas )

#### Specifying Model Using C Snippets ####

