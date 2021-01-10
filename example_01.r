loc <- url("https://kingaa.github.io/sbied/intro/parus.csv")
dat <- read.csv(loc)
head(dat)

plot(pop~year,data=dat,type='o')

library(pomp)
parus <- pomp(dat,times="year",t0=1959)
plot(parus)


# Dynamics 
stochStep <- Csnippet("
  e = rnorm(0,sigma);
  N = r*N*exp(-c*N+e);
")

# Adds dynamics to pomp defined instance
# Define time step delta.t to 1 year
pomp(parus,rprocess=discrete_time(step.fun=stochStep,delta.t=1),
     paramnames=c("r","c","sigma"),
     statenames=c("N","e")) -> parus

sim <- simulate(parus,
                params=c(N.0=1,e.0=0,r=12,c=1,sigma=0.5),
                format="data.frame")

plot(N~year,data=sim,type='o')

# Measurement model
rmeas <- Csnippet("pop = rpois(phi*N);")

# Likelihood
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")

# Add measurement Model and Likelihood to pomp object
pomp(parus,
     rmeasure=rmeas,
     dmeasure=dmeas,
     statenames=c("N"),
     paramnames=c("phi")) -> parus

coef(parus) <- c(N.0=1, e.0=0, r=20, c=1, sigma=0.1,phi=200)

library(ggplot2)
sims <- simulate(parus,nsim=3,format="d",include.data=TRUE)
ggplot(data=sims,mapping=aes(x=year,y=pop))+geom_line()+
  facet_wrap(~.id,ncol=1,scales="free_y")

skel <- Csnippet("DN = r*N*exp(-c*N);")
parus <- pomp(parus,skeleton=map(skel),paramnames=c("r","c"),statenames=c("N"))

traj <- trajectory(parus,params=c(N.0=1,r=12,c=1),format="data.frame")
plot(N~year,data=sim,type='o')
lines(N~year,data=traj,type='l',col='red')


#### Simulate the model ####
x <- simulate(parus)
plot(x)

#### Multiple Simulations ####
x <- simulate(parus,nsim=9,format="data.frame",include.data=TRUE)
ggplot(data=x,aes(x=year,y=pop,group=.id,color=(.id=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(~.id,ncol=2)

#### Change parameters ####
coef(parus,c("r","N.0","sigma")) <- c(2,1,0.5)
coef(parus)
x = simulate(parus,format="data.frame",include.data=TRUE)
ggplot(data=x,aes(x=year,y=pop,group=.id,color=(.id=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(1,ncol=1)

