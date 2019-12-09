## #######################################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    December 2019
## Notes:   Purpose is to replicate the SD model from Alonge 2017
##          Paper here: https://academic.oup.com/heapol/article/32/10/1417/4210153
##          Secondarily, will perfrom new sensitivity analyses
## #######################################################################


# load libraries
libs <- c('deSolve','data.table','ggplot2')
for(l in libs) library(l,character.only = TRUE)

#  simulation times and time step and time vector
START   <- 0
FINISH  <- 200
STEP    <- 1
simtime <- seq(START, FINISH, by=STEP)

# set up initial stock values for Revenue, Quality, Volume of Services
stocks <- c(R = 0, Q  = 0, V = -0.5)

# model parameters (baseline scenario)
auxs   <- c(rR        = 0.3, 
            rRQ       = 0.2,
            rQ        = 0.3, 
            rV        = 0.1, 
            Qdelay    = 3,
            Vdelay    = 0,
            P4Pfactor = 0,
            G_0       = 0,
            M_0       = 0)

# Logistic forcing function used in paper to limit stocks from -1 to 1
L <- function(x) { (1-exp(-2*x))/(1+exp(-2*x))  }

# The Model function, takes 3 arguments from ode()
model <- function(time, stocks, auxs){
  with(as.list(c(stocks, auxs)),{
    
    # delay variables
    if (time >= Qdelay)
      Qlag <- lagvalue(time - Qdelay, 2)
    else
      Qlag <- lagvalue(1,2)
    
    if (time > Vdelay)
      Vlag <- lagvalue(time - Vdelay, 3)
    else
      Vlag <- lagvalue(1,2)
    
    # dependent variables
    P4P <- P4Pfactor * ifelse(L(Vlag)<0,0,L(Vlag))
    G   <- G_0 + 1.1^P4P - 1
    M   <- M_0 + P4P - G
    
    # differential equations
    dR_dt <- rR*L(V) - rRQ*R*(1.3^M) + P4P
    
    if(V>0)
      dQ_dt <- rRQ*R*(1.3^M) - rQ*L(V)*(1.1^G)/(1.3^M)
    else 
      dQ_dt <- rRQ*R*(1.3^M) - rQ*L(V)*(1.3^M)/(1.1^G)
    
    if(Qlag>0)
      dV_dt <- rV*L(Qlag)*(1.3^M)
    else 
      dV_dt <- rV*L(Qlag)/(1.3^M)
    
    # All the results for the time step
    ans <- list(c(dR_dt,dQ_dt,dV_dt))
  })
}

## replicate figures 5 and 6 in alonge et al paper by changing the aux variables

# helper function to replace parameter values from baseline
prepl <- function(x, p = auxs){
  if(!is.null(x))
    for(newp in names(x)){
      p[[newp]] <- x[[newp]]
    }
  return(p)
}

# replace values according to Table A3 from the appendix:
aux1  <- prepl(x = NULL)
aux2  <- prepl(c(P4Pfactor = 0.2))
aux3a <- prepl(c(P4Pfactor = 0.2, M_0 = 0.3))
aux3b <- prepl(c(P4Pfactor = 0.2, G_0 = 0.2))
aux3c <- prepl(c(P4Pfactor = 0.2, G_0 = 0.8))
aux3d <- prepl(c(P4Pfactor = 0.2, Vdelay = 3))
aux3e <- prepl(c(P4Pfactor = 0.081))
aux4a <- prepl(c(P4Pfactor = 0.2, Vdelay = 3, M_0 = 0.1,  G_0 = 0.3))
aux4b <- prepl(c(P4Pfactor = 0.2, Vdelay = 3, M_0 = 0.05, G_0 = 0.8))
aux4c <- prepl(c(P4Pfactor = 0.2, Vdelay = 3, M_0 = 0.3,  G_0 = 0.05))

# solve models and plot for fig 5:
o5 <- do.call('rbind', list(
      melt(data.table(dede(y=stocks, times=simtime, func = model, parms=aux1)), id='time'),
      melt(data.table(dede(y=stocks, times=simtime, func = model, parms=aux2)), id='time'),
      melt(data.table(dede(y=stocks, times=simtime, func = model, parms=aux3a)), id='time'),
      melt(data.table(dede(y=stocks, times=simtime, func = model, parms=aux3c)), id='time') ))
o5$label <- rep(c('a) Baseline','b) P4P Only','c) P4P w/ Motivation','d) P4P w/ High Gaming'), 
               each = nrow(o5)/4)

ggplot(o5[variable!='R'], aes(x=time,y=L(value),color=variable)) + geom_line(lwd=1.2) +
  ylim(-1,1) + xlim(0,100) + theme_bw() + geom_hline(yintercept = 0, lty = 'dotted') +
  scale_color_manual(values=c("#2EC4B6", "#E71D36"), labels = c('Quality', "Volume")) +
  facet_wrap(~label) + ggtitle('Alonge et al. 2017 Fig. 5 Replication') +
  theme(legend.position = 'bottom') + labs(color = '')


# solve models and plot for fig 6:
o6 <- do.call('rbind', list(
  melt(data.table(dede(y=stocks, times=simtime, func = model, parms=aux4a)), id='time'),
  melt(data.table(dede(y=stocks, times=simtime, func = model, parms=aux4b)), id='time'),
  melt(data.table(dede(y=stocks, times=simtime, func = model, parms=aux4c)), id='time') ))
o6$label <- rep(c('a) Equal allocation','b) Proportionate to salaries','c) Proportionate to services'), 
               each = nrow(o6)/3)

ggplot(o6[variable!='R'], aes(x=time,y=L(value),color=variable)) + geom_line(lwd=1.2) +
  ylim(-1,1) + xlim(0,100) + theme_bw() + geom_hline(yintercept = 0, lty = 'dotted') +
  scale_color_manual(values=c("#2EC4B6", "#E71D36"), labels = c('Quality', "Volume")) +
  facet_wrap(~label,ncol=2) + ggtitle('Alonge et al. 2017 Fig. 5 Replication') +
  theme(legend.position = 'bottom') + labs(color = '')



