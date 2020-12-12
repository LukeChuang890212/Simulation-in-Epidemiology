rm(list=ls())

setwd("C:/Users/莊明儒/Desktop/Epidemiology")
source("Functions.r")

N = 50:200 #no. of samples 

p.mr = 0.02  #mortality rate in pupulation 
pers.m.in.ltfus = seq(0.0,0.1,0.02) #percentages of mortality in loss-to-follow-up subjects (relatedness of loss-to-follow-up to mortality)

n.sim = 1000 #no. of simulation (fixed)
ltfurs = seq(0.0,0.5,0.1) #rate of loss to follow up (fixed)
n.followyr = 10 #no. of years to follow up (fixed)

n.N = length(N) # no. of differnt sample sizes

for(ltfur in ltfurs){
  avg.ltfur = 1-(1-ltfur)^(1/n.followyr)
  
  reses = list()
  
  png(paste("emr~ss(",ltfur,',',min(pers.m.in.ltfus),'~',max(pers.m.in.ltfus),').png',sep=''),width=1436,height=581)
  par(mfrow=c(2,3))
  for(per.m.in.ltfus in pers.m.in.ltfus){
    res = list(m.mrs=rep(0,n.N),downlm.ci=rep(0,n.N),uplm.ci=rep(0,n.N),bias=rep(0,n.N))
    for(i in 1:n.N){ # i means index of sample size 
      n=N[i] # sample size
        
      mrs = rep(0,n.sim)
      # n.total.ltfuss = rep(0,n.sim)
      for(sim in 1:n.sim){ # sim means the (sim)th time of simulation
        subjects = sample(n)
        mrs[sim] = cal.mr(subjects,p.mr,per.m.in.ltfus,n.followyr,avg.ltfur) 
      }
    
      res = cal.m.ci(i,res,mrs) # calculate the mean and ci for our estimation of mortality rate and store them into res list at index i
      res = cal.bias(i,p.mr,res,mrs) # calculate the bias for our estimation of mortality rate and store them into res list at index i
    }
    reses[[length(reses)+1]] = res
    print(paste("Finished simulation for level of association =",per.m.in.ltfus,"\nunder proportion of loss of follow-up =",ltfur))
    
    #plot(x=N,y=res$bias,type='l',xlab="Sample Size",ylab="Bias(rmse)",main=paste("Level of Association =",per.m.in.ltfus)) # plot bias~ss
    plot.emr.ss(p.mr,per.m.in.ltfus,res) # plot emr~ss
  }
  #title(paste("Proportion of Loss of Follow-up =",ltfur),outer=T)
  
  dev.off()
  
  png(paste("bias~ss(",ltfur,',',min(pers.m.in.ltfus),'~',max(pers.m.in.ltfus),').png',sep=''),width=1436,height=581)
  par(mfrow=c(2,3))
  for(i in 1:length(pers.m.in.ltfus)){ # i means index of per.m.in.ltfus
    per.m.in.ltfus = pers.m.in.ltfus[i]
    plot(x=N,y=reses[[i]]$bias,type='l',xlab="Sample Size",ylab="Bias(rmse)",main=paste("Level of Association =",per.m.in.ltfus)) # plot bias~ss
  }
  #title(paste("Proportion of Loss of Follow-up =",ltfur),outer=T)
  
  dev.off()
}












# test for simulating mean of mr 
cal.mr = function(subjects,p.mr,per.m.in.ltfus,n.followyr,avg.ltfur){
  n.total.mos = 0 # no. of total mortality case
  # n.total.ltfus = 0
  ptar = 0 #population time at risk
  for(yr in 1:n.followyr){
    n.mos = rpois(n=1,lambda=length(subjects)*p.mr) # no. of mortality-occurence subjects
    n.ltfus = rbinom(n=1,size=length(subjects)+n.total.mos,prob=avg.ltfur) # no. of ltfus
    # n.total.ltfus = n.total.ltfus+n.ltfus
    # n.total.ltfuss[sim] = n.total.ltfus
    # 
    # print(paste("n.ltfus:",n.ltfus))
    # print(paste("n.mos:",n.mos))
    
    ltfus = sample(subjects,n.ltfus) # loss to follow-up subjects
    
    subjects = subjects[! subjects %in% ltfus] # remain only subjects who is not ltfus (1st time to del subjects)
    
    nominal.n.mos = n.mos-round(n.ltfus*per.m.in.ltfus) # n.ltfus*per.m.in.ltfus = no. of subjects who is ltfus due to mortality)
    if(nominal.n.mos >0){ # if there is any extra dropout due to mortality
      nominal.mos = sample(subjects,nominal.n.mos) # nomial no. of mos 
      
      subjects = subjects[! subjects %in% nominal.mos] # remain only subjects who is not nominal.mos (2nd time to del subjects)
      
      n.total.mos = n.total.mos+nominal.n.mos
    }
    
    ptar = ptar+n.ltfus*(yr)+nominal.n.mos*(yr-0.5)   # add the contribution of person-year from each subject in ltfus & mos to ptar
  }
  ptar = ptar+length(subjects)*n.followyr # add the contribution of person-year from each remain subjects to ptar
  
  # print(paste("length(subjects):",length(subjects)))
  # print(paste("ptar:",ptar))
  # print(paste("n.total.mos:",n.total.mos))
  # print(paste("n.total.ltfus:",n.total.ltfus))
  
  mr = n.total.mos/ptar
  # print(paste("mr:",mr))
  
  return(mr)
}

uplm.ci = round(mean(mrs)+1.96*sd(mrs),4) # uplimit of ci
downlm.ci = round(mean(mrs)-1.96*sd(mrs),4) # downlimit of ci
downlm.ci = ifelse(downlm.ci<0,0,downlm.ci) # if downlimit is negative then replace it with 0 
print(paste("CI: [",downlm.ci,",",uplm.ci,"]"))

# test for avg.ltfur
n=200
subjects = sample(n)
losses = c()
for(j in 1:1000){
    len = length(subjects)
  loss = 0
  for(i in 1:10){
    loss = loss+rbinom(n=1,size=len,prob=avg.ltfur)
    len = len-rbinom(n=1,size=len,prob=avg.ltfur)
  }
  print(loss)
  losses = c(losses,loss)
}
mean(losses)

#Note
#54th line