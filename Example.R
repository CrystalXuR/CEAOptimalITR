#--------------------------------------------------------------------------------------------------------------#######
#-------------------- This is an example of estimating the most cost-effective         ------------------------#######  
#-------------------- individualized treatment rule using the "OptimalCEITR" function  ------------------------#######
#--------------------------------------------------------------------------------------------------------------#######
#setwd("PATH of the saved OptimalCEITR function and example data")
source("./OptimalCEITR.R")

# Load example data 
dat <- read.csv("./ExampleData.csv")

# Estimated censoring weights
dat$complete <- ifelse((dat$event==1 | dat$FT>=dat$tau),1,0)   #complete:having event|contribute full follow-up of interest
cen.wm       <- glm(complete ~ x1+x2+x3+A+A:x1+A:x2+A:x3, family = "binomial", data=dat)
dat$cen.w    <- predict(cen.wm, type = "response")

# Estimated individualized treatment rule 
Xs <- dat[,1:3]
OptimalITR <- OptimalCEITR(X=Xs,Xs=Xs,Ms=dat$x2,A=dat$A,ST=dat$FT,event=dat$event,C=dat$CC,cen=dat$complete,
                           est_cen=dat$cen.w,tau=dat$tau[1],lambda=dat$lambda[1],data=dat)
Est.g.opt <- as.factor(OptimalITR[[1]])
VIP <- OptimalITR[[2]] 

# Assess performance
datcomp <- dat[dat$complete==1,]
datcomp$g.opt <- as.factor(datcomp$g.opt)
table(datcomp$g.opt,Est.g.opt)    #CCR = (158+597)/dim(datcomp)[1]=85.4%

library(ggplot2);library(gridExtra)
g.opt.DB <- qplot(x1, x2, colour = g.opt, data=datcomp)+
            scale_colour_manual(values=c("#0000CC", "#FF3300"))+
            ggtitle("g.opt")+
            theme(plot.title = element_text(size=12),
                  axis.text = element_text(size=12),axis.title = element_text(size=12))

Est.g.opt.BD <- qplot(x1, x2, colour = Est.g.opt, data=datcomp)+
                scale_colour_manual(values=c("#0000CC", "#FF3300"))+
                ggtitle("Est.g.opt")+
                theme(plot.title = element_text(size=12),
                axis.text = element_text(size=12),axis.title = element_text(size=12))

print(grid.arrange(arrangeGrob(g.opt.DB, Est.g.opt.BD, nrow=1, ncol=2), nrow=1, heights=c(10)))  #Decision Boundaries
##-- End of the Code
