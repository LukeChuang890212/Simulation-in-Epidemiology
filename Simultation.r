setwd("C:/Users/莊明儒/Desktop/Epidemiology")

N = 50:200 #no. of samples 

p.mr = 0.02  #mortality rate in pupulation 
per.m.in.ltfus = 0.0 #percentage of mortality in loss-to-follow-up subjects (relatedness of loss-to-follow-up to mortality)

n.sim = 1000 #no. of simulation (fixed)
ltfur = 0.5 #rate of loss to follow up (fixed)
n.followyr = 10 #no. of years to follow up (fixed)
avg.ltfur = 1-(1-ltfur)^(1/n.followyr)

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

n=500
mrs = rep(0,n.sim)
# n.total.ltfuss = rep(0,n.sim)
for(sim in 1:n.sim){
  subjects = sample(n)
  mrs[sim] = cal.mr(subjects,p.mr,per.m.in.ltfus,n.followyr,avg.ltfur) 
}
mean(mrs)

uplm.ci = round(mean(mrs)+1.96*sd(mrs),4) # uplimit of ci
downlm.ci = round(mean(mrs)-1.96*sd(mrs),4) # downlimit of ci
downlm.ci = ifelse(downlm.ci<0,0,downlm.ci) # if downlimit is negative then replace it with 0 
print(paste("CI: [",downlm.ci,",",uplm.ci,"]"))

# test for simulate mean of mr 
n=500
mrs = rep(0,n.sim)
# n.total.ltfuss = rep(0,n.sim)
for(sim in 1:n.sim){
  subjects = sample(n)
  
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
  
  mrs[sim] = mr 
}
mean(mrs)
# (n.total.ltfuss)

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