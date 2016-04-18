data_one <- read.table("data_onesided.txt")
data_two <- read.table("data_twosided.txt")
data_neut <- read.table("data_neutral.txt")


data_onesided <- matrix(NA,nrow=nrow(data_one),ncol=3)
data_onesided[,1] <- data_one[,1]

data_onesided_two <- matrix(NA,nrow=nrow(data_one),ncol=3)
data_onesided_two[,1] <- data_one[,1]

data_twosided <- matrix(NA,nrow=nrow(data_two),ncol=3)
data_twosided[,1] <- data_two[,1]

data_neutral <- matrix(NA,nrow=nrow(data_neut),ncol=3)
data_neutral[,1] <- data_neut[,1]

# for onesided vs twosided
data_onesided_two[,2] <- rep(1,nrow(data_onesided_two))
data_onesided_two[,3] <- rep(0,nrow(data_onesided_two))

data_twosided[,2] <- rep(0,nrow(data_twosided))
data_twosided[,3] <- rep(1,nrow(data_twosided))

data_orig_one_two <- rbind(data_onesided_two,data_twosided)

colnames(data_orig_one_two) <- c("ambivalence","onesided","twosided")
write.table(data_orig_one_two, sep="\t", file="Data_orig_onesided_twosided.txt",row.names = F)

t_one_two <- t.test(data_onesided_two[,1],data_twosided[,1],var.equal=T)

# for neutral vs onesided
data_neutral[,2] <- rep(1,nrow(data_neutral))
data_neutral[,3] <- rep(0,nrow(data_neutral))

data_onesided[,2] <- rep(0,nrow(data_onesided))
data_onesided[,3] <- rep(1,nrow(data_onesided))

data_orig_neut_one <- rbind(data_neutral, data_onesided)

colnames(data_orig_neut_one) <- c("ambivalence","neutral","onesided")
write.table(data_orig_neut_one, sep="\t", file="Data_orig_neutral_onesided.txt",row.names = F)

t_neut_one <- t.test(data_neutral[,1],data_onesided[,1],var.equal=T)

# for neutral vs twosided
data_orig_neut_two <- rbind(data_neutral,data_twosided)

colnames(data_orig_neut_two) <- c("ambivalence","neutral","twosided")
write.table(data_orig_neut_two, sep="\t", file="Data_orig_neutral_twosided.txt",row.names = F)

t_neut_two <- t.test(data_neutral[,1],data_twosided[,1],var.equal=T)

