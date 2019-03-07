setwd('C:/Users/Asus/Documents/p04_prs/')
setwd('p01_examine_sumstats/')
sum_stats = read.table("C:/Users/Asus/Documents/p04_prs/data/DIAGRAMv3.2012DEC17.txt", header = T, as.is=T, nrows=5000)
head(sum_stats)
names(sum_stats)<-c("SNP","chr","pos","a1","a0","p","or1","or1_95l","or1_95u","n_cases","n_controls")    
p_thresh_for_all<-1e-3
sum_stats0 <- sum_stats[sum_stats$p < p_thresh_for_all,] 
manhattan(sum_stats0[sum_stats0$p>1e-20,],chr="chr",p="p",bp="pos")
