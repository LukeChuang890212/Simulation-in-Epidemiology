cal.mr = function(subjects,p.mr,per.m.in.ltfus,n.followyr,avg.ltfur){
  n.total.mos = 0 # no. of total mortality case
  # n.total.ltfus = 0
  ptar = 0 #population time at risk
  for(yr in 1:n.followyr){
    # print(paste("start length(subjects):",length(subjects)))
    n.mos = rpois(n=1,lambda=length(subjects)*p.mr) # no. of mortality-occurence subjects
    n.ltfus = rbinom(n=1,size=length(subjects)+n.total.mos,prob=avg.ltfur) # no. of ltfus
    # n.total.ltfus = n.total.ltfus+n.ltfus
    # n.total.ltfuss[sim] = n.total.ltfus
    
    # print(paste("n.ltfus:",n.ltfus))
    # print(paste("n.mos:",n.mos))
    if(n.ltfus<=length(subjects)){ # if no. of loss of follow-up is still smaller than th remain subjects
      ltfus = sample(subjects,n.ltfus) # loss to follow-up subjects
    }  
    
    subjects = subjects[! subjects %in% ltfus] # remain only subjects who is not ltfus (1st time to del subjects)
    
    
    nominal.n.mos = round(n.mos-n.mos*per.m.in.ltfus) # n.mos*per.m.in.ltfus = no. of subjects who is ltfus due to mortality) # where rounding error might happen
    # print(paste("nominal.n.mos:",nominal.n.mos))
    
    if(nominal.n.mos<=length(subjects)){ # if no. of extra mortality cases is still smaller than th remain subjects
      nominal.mos = sample(subjects,nominal.n.mos) # nomial no. of mos 
      
      subjects = subjects[! subjects %in% nominal.mos] # remain only subjects who is not nominal.mos (2nd time to del subjects)
      
      n.total.mos = n.total.mos+nominal.n.mos
    }
    
    ptar = ptar+n.ltfus*(yr)+nominal.n.mos*(yr-0.5)   # add the contribution of person-year from each subject in ltfus & mos to ptar
  }
  ptar = ptar+length(subjects)*n.followyr # add the contribution of person-year from each remain subjects to ptar
  
  # print(paste("end length(subjects):",length(subjects)))
  # print(paste("ptar:",ptar))
  # print(paste("n.total.mos:",n.total.mos))
  # print(paste("n.total.ltfus:",n.total.ltfus))
  
  mr = n.total.mos/ptar
  # print(paste("mr:",mr))
  
  return(mr)
}

cal.m.ci = function(i,res,mrs){
  m.mrs = mean(mrs)
  sd.mrs = sd(mrs)
  print(paste("mean:",m.mrs))
  
  uplm.ci = round(m.mrs+1.96*sd.mrs,4) # uplimit of ci
  downlm.ci = round(m.mrs-1.96*sd.mrs,4) # downlimit of ci
  downlm.ci = ifelse(downlm.ci<0,0,downlm.ci) # if downlimit is negative then replace it with 0 
  print(paste("CI: [",downlm.ci,",",uplm.ci,"]"))
  
  res$m.mrs[i] = m.mrs
  res$downlm.ci[i] = downlm.ci
  res$uplm.ci[i] = uplm.ci
  
  return(res)
}

cal.bias = function(i,p.mr,res,mrs){
  #bias = sum(abs(mrs-p.mr)/p.mr)/length(mrs)
  bias = sum(mrs-p.mr)/length(mrs)
  # bias = sqrt(sum((mrs-p.mr)^2)/length(mrs))
  res$bias[i] = bias
  
  return(res)
}

plot.emr.ss = function(p.mr,per.m.in.ltfus,res){
  bp = barplot(height=res$m.mrs,names.arg=N,xlab="Sample Size",ylab="Estimated Mortality Rate",ylim=c(0.0,max(res$uplm.ci)),main=paste("Level of Association =",per.m.in.ltfus))
  lines(bp,res$downlm.ci,lty=2)
  lines(bp,res$uplm.ci,lty=2)
  abline(h=p.mr,lwd=2)
}

