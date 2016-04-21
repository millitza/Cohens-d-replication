#### THIS DATA FILE WAS BASED ON THE DATA FILE PROVIDED BY THE OPEN SCIENCE FRAMEWORK. ####
##@@ GENERAL INFORMATION @@##
#@ Study Title: Replication of The Effects of an Implemental Mind-set on Attitude Strength by Henderson, de Liver, & Gollwitzer (2008, Journal of Personality and Social Psychology) 
#@ Coder name: Julian Mestwerdt
#@ Coder e-mail: julian.mestwerdt@gmail.com
#@ Type of statistic: F
#@ Type of effect-size: Partial eta squared
#@ OSF link project: https://osf.io/79dey/
#@ OSF link replication report: https://osf.io/cjr7d/

##@@ DATA LOADING @@##
#@ NOTE: Data must be loaded from OSF directly and NOT rely on any local files.
Data <- read.delim("reconciling.for.repro.txt", header=TRUE)
##@@ DATA MANIPULATION @@##
#@ NOTE: Include here ALL difference between OSF data and data used in analysis
#@ TIP: You will want to learn all about dplyr for manipulating data.
Data <- Data[Data$Include==1,] # Excluding cases
Data$Indecisive <- as.numeric(as.character(Data$Indecisive))
Data$Condition  <- as.factor(Data$Condition)

Data$NotTorn     <- Data$NotTorn*-1   # Reverse code item 2
Data$Indecisive  <- Data$Indecisive-4 # Equating the scoring of item 3 (1 to 7) to 
# the scoring of item 1 & 2 (-3 to 3) 

Data$value <- (Data$StrongMixed + Data$NotTorn + Data$Indecisive)/3
# Creating the ambivalence score, based on the mean of item 1 to 3                                                             

Data[49,33] <- (Data[49,22]+Data[49,22])/2
# This participant did not provide an answer for item 3.
# The ambivalence score for this participant is calculated to be the mean of item 1 & 2.

Data_essential <- cbind(Data$value,Data$Condition)
colnames(Data_essential) <- c("ambivalence","condition")

Ncond <- table(Data_essential[,2])

#one-sided vs two-sided
Data_exp_one_two <- matrix(NA,nrow=Ncond[2]+Ncond[3],ncol=3)
Data_exp_one_two[,1] <- Data_essential[25:70,1]
Data_exp_one_two[,2] <- c(rep(1,Ncond[2]),rep(0,Ncond[3]))
Data_exp_one_two[,3] <- c(rep(0,Ncond[2]),rep(1,Ncond[3]))
write.table(Data_exp_one_two, sep="\t",file="Data_repl_onesided_twosided.txt",row.names=F)

t_onesided_twosided <- t.test(Data_exp_one_two[1:23,1],Data_exp_one_two[24:46,1],var.equal=T)

#neutral vs one-sided
Data_exp_neut_one <- matrix(NA,nrow=Ncond[1]+Ncond[2],ncol=3)
Data_exp_neut_one[,1] <- c(Data_essential[1:24,1],Data_essential[25:47,1])
Data_exp_neut_one[,2] <- c(rep(1,Ncond[1]),rep(0,Ncond[2]))
Data_exp_neut_one[,3] <- c(rep(0,Ncond[1]),rep(1,Ncond[2]))
write.table(Data_exp_neut_one, sep="\t",file="Data_repl_neutral_onesided.txt",row.names=F)

t_neutral_onesided <- t.test(Data_exp_neut_one[1:24,1],Data_exp_neut_one[25:47,1],var.equal=T)

#neutral vs two-sided
Data_exp_neut_two <- matrix(NA,nrow=Ncond[1]+Ncond[3],ncol=3)
Data_exp_neut_two[,1] <- c(Data_essential[1:24,1],Data_essential[48:70,1])
Data_exp_neut_two[,2] <- c(rep(1,Ncond[1]),rep(0,Ncond[3]))
Data_exp_neut_two[,3] <- c(rep(0,Ncond[1]),rep(1,Ncond[3]))
write.table(Data_exp_neut_two, sep="\t",file="Data_repl_neutral_twosided.txt",row.names=F)

t_neutral_twosided <- t.test(Data_exp_neut_two[1:24,1],Data_exp_neut_two[25:47,1],var.equal=T)


