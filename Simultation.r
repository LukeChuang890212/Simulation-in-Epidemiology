rm(list=ls())

setwd("C:/Users/莊明儒/Desktop/Epidemiology")
source("Functions.r")

N = 50:200 #no. of samples 

p.mrs = c(0.11) #mortality rate in pupulation c(0.02,0.1,0.4) 
pers.m.in.ltfus = c(0.1) #percentages of mortality that is mistaken for loss of follow-up (association of loss-to-follow-up with mortality)

n.sim = 1000 #no. of simulation (fixed)
ltfurs = seq(0.1,0.5,0.1) #rate of loss to follow up (fixed)
n.followyr = 10 #no. of years to follow up (fixed)

n.N = length(N) # no. of differnt sample sizes

for(p.mr in p.mrs){
  save.dir = paste("bias.x.emr~ss(",p.mr,')',sep='')
  dir.create(save.dir,showWarnings = FALSE)
  
  # biases = list() # store the biases for each ltfur (proportion of loss of follow up)
  for(ltfur in ltfurs){
    avg.ltfur = 1-(1-ltfur)^(1/n.followyr)
    
    reses = list() #save each res of different per.m.in.lfus for the plotting of bias~ss into reses list 
    
    png(paste(save.dir,'/',"mr~ss(",ltfur,',',min(pers.m.in.ltfus),'~',max(pers.m.in.ltfus),').png',sep=''),width=957,height=581)
    par(mfrow=c(2,2))
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
        
        #hist(mrs,breaks="Scott")
        #Sys.sleep(60)
        
        res = cal.m.ci(i,res,mrs) # calculate the mean and ci for our estimation of mortality rate and store them into res list at index i
        res = cal.bias(i,p.mr,res,mrs) # calculate the bias for our estimation of mortality rate and store them into res list at index i
      }
      reses[[length(reses)+1]] = res
      cat(paste("Simulation for level of association =",per.m.in.ltfus,"under \nproportion of loss of follow-up =",ltfur,'&','population inccidence rate =',p.mr,"\nhas been done"))
      
      # plot(x=N,y=res$bias,type='l',xlab="Sample Size",ylab="Bias(rmse)",main=paste("Level of Association =",per.m.in.ltfus)) # plot bias~ss
      plot.emr.ss(p.mr,per.m.in.ltfus,res) # plot emr~ss
    }
    #title(paste("Proportion of Loss of Follow-up =",ltfur),outer=T)
    
    dev.off()
    
    png(paste(save.dir,'/',"bias~ss(",ltfur,',',min(pers.m.in.ltfus),'~',max(pers.m.in.ltfus),').png',sep=''),width=957,height=581)
    par(mfrow=c(2,2))
    for(i in 1:length(pers.m.in.ltfus)){ # i means index of per.m.in.ltfus
      per.m.in.ltfus = pers.m.in.ltfus[i]
      plot(x=N,y=reses[[i]]$bias,type='l',xlab="Sample Size",ylab="Bias",main=paste("Level of Association =",per.m.in.ltfus)) # plot bias~ss
    }
    #title(paste("Proportion of Loss of Follow-up =",ltfur),outer=T)
    
    dev.off()
    
    # biases[[length(biases)+1]] = reses[[length(reses)]]$bias
  }
  
  # png(paste(save.dir,'/',"bias~ss(proportion=",min(ltfurs),'~',max(ltfurs),').png',sep=''),width=957,height=581)
  # par = (mfrow=c(1,1))
  # lcols = c("darkred","darkorange","darkgreen","darkblue","purple") # lines colors
  # for(ltfur in ltfurs){
  #   bias.plot.index = match(ltfur,ltfurs)
  #   if(bias.plot.index == 1){
  #     plot(x=N,y=biases[[bias.plot.index]],type='l',ylim=c(min(sapply(biases, min)),max(sapply(biases, max))),xlab="Sample Size",ylab="Bias",main=paste("Level of Association =",per.m.in.ltfus),col=lcols[bias.plot.index]) # plot bias~ss
  #   }else{
  #     lines(x=N,y=biases[[bias.plot.index]],lty=1,col=lcols[bias.plot.index])
  #     if(bias.plot.index == length(ltfurs)){
  #       legend("topright",legend=c(paste("proportion =",ltfurs)),lty=rep(1,length(ltfurs)),col=lcols)
  #     }
  #   }
  # }
  # 
  # dev.off()
}

# 隨著母體參數越大，低估效果隨樣本數增加而變大的現象越明顯。
# 低估效果只在關聯程度為0.08或0.1的時候才有。隨關聯程度越高低估效果亦越明顯。










