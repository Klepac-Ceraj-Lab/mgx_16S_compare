library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(sna)
library(reshape2)
library(gtools)
library(ape)

setwd("/Users/danielle/Documents/thesis/paper-method-comparison")

df <- read.csv("transposed_mgxamp_df.csv", header=TRUE)

total_columns <- ncol(df)
total_samples <- nrow(df)

abund_table <- as.matrix(df[,5:total_columns])
bins = seq(min(abund_table), max(abund_table), by=0.01)
hist(abund_table[abund_table != 0.0], xlim = c(0,0.01), col = 'skyblue3',
     breaks = 10000)

# calculating Shannon diversity (alpha-diversity)
df["shannon"] <- diversity(abund_table, "shannon")


# adding developmental stage

df$dev_stage[df$AgeMonths<15 ] <- "less than 15 months"
df$dev_stage[(df$AgeMonths>=15 & df$AgeMonths <30) ] <- "15 to 30 months"
df$dev_stage[df$AgeMonths>30 ] <- "older than 30 months"

# finding number of kids in each developmental stage
length(df$dev_stage[df$AgeMonths<15])/2
length(df$dev_stage[(df$AgeMonths>=15 & df$AgeMonths <30) ])/2
length(df$dev_stage[df$AgeMonths>30])/2

df$dev_stage <- factor(df$dev_stage,
                       levels = c("less than 15 months", 
                                  "15 to 30 months", 
                                  "older than 30 months"), ordered = TRUE)

# Plot 1: plot shannon diversity, broken down by age and profiling method

p1 <- ggplot(df, aes(x = dev_stage, y = shannon)) + 
  geom_boxplot(aes(colour = method))
p1 <- p1 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
  ylab("Shannon diversity")+labs(title = "", tag = "A") + 
  xlab("developmental stage") +
  theme(legend.position = "none")

p1

# shannon diversity for different ages based on profiling method 
# 16S

df_amp <- df[df$method=="amp",]
mean(df_amp[df_amp$dev_stage == "less than 15 months",]$shannon, na.rm=TRUE)
mean(df_amp[df_amp$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
mean(df_amp[df_amp$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)

# comparing shannon index by developmental stage
t.test((df_amp[df_amp$dev_stage == "less than 15 months",]$shannon),
       df_amp[df_amp$dev_stage == "15 to 30 months",]$shannon)
t.test((df_amp[df_amp$dev_stage == "15 to 30 months",]$shannon),
       df_amp[df_amp$dev_stage == "older than 30 months",]$shannon)
t.test((df_amp[df_amp$dev_stage == "less than 15 months",]$shannon),
       df_amp[df_amp$dev_stage == "older than 30 months",]$shannon)
# comparing kids older and younger than 30 months 
t.test((df_amp[df_amp$dev_stage == "less than 15 months" | df_amp$dev_stage =="15 to 30 months",]$shannon),
       df_amp[df_amp$dev_stage == "older than 30 months",]$shannon)


# mgx
df_mgx <- df[df$method=="mgx",]
mean(df_mgx[df_mgx$dev_stage == "less than 15 months",]$shannon, na.rm=TRUE)
mean(df_mgx[df_mgx$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
mean(df_mgx[df_mgx$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)

# comparing shannon index by developmental stage
t.test((df_mgx[df_mgx$dev_stage == "less than 15 months",]$shannon),
       df_mgx[df_mgx$dev_stage == "15 to 30 months",]$shannon)
t.test((df_mgx[df_mgx$dev_stage == "15 to 30 months",]$shannon),
       df_mgx[df_mgx$dev_stage == "older than 30 months",]$shannon)
t.test((df_mgx[df_mgx$dev_stage == "less than 15 months",]$shannon),
       df_mgx[df_mgx$dev_stage == "older than 30 months",]$shannon)

# comparing kids older and younger than 30 months 
t.test((df_mgx[df_mgx$dev_stage == "less than 15 months" | df_mgx$dev_stage =="15 to 30 months",]$shannon),
       df_mgx[df_mgx$dev_stage == "older than 30 months",]$shannon)


anova <- aov(shannon ~ dev_stage, data = df)
TukeyHSD(anova)

# paired t-test
t.test((df[df$dev_stage == "less than 15 months" & 
             df$method == "amp",]$shannon),
       df[df$dev_stage == "less than 15 months" & 
            df$method == "mgx",]$shannon,paired=TRUE)
t.test((df[df$dev_stage == "15 to 30 months" & 
             df$method == "amp",]$shannon),
       df[df$dev_stage == "15 to 30 months" & 
            df$method == "mgx",]$shannon,paired=TRUE)
t.test((df[df$dev_stage == "older than 30 months" & 
             df$method == "amp",]$shannon),
       df[df$dev_stage == "older than 30 months" & 
            df$method == "mgx",]$shannon,paired=TRUE) ## this number is in the paper


# bray curtis dissimilarity for samples

# BC distance for kids of different ages
### bug_start = first column of data frame that includes relative abundance data
### bug_end = last column of data frame that incldues relative abundance data
### dev_stage_arg = developmental stage (<15, 15-30, >30 months)
### method_arg = method to profile microbial community (16S or metagenomics)
### returns a data frame of bray-curtis dissimilarity values for each age 
### group/profiling method

calc_bc_dist <- function(bug_start, bug_end, dev_stage_arg, method_arg) {
  df_slice <- df[df$dev_stage == dev_stage_arg
                 & df$method == method_arg,]
  
  abund_table <- df_slice[,bug_start:bug_end]
  abund_table <- abund_table[,colSums(abund_table) > 0]
  bc_matrix <- as.matrix(vegdist(abund_table, "bray", diag=FALSE,upper=FALSE))
  bc_matrix[bc_matrix == 0.0] <- NA
  bc_vec <- na.omit(as.vector(bc_matrix))
  bc_df <- data.frame(bc_vec)
  colnames(bc_df)<-c("BC_dist")
  bc_df["method"] <- method_arg
  bc_df["dev_stage"] <- dev_stage_arg
  
  return(bc_df)
}

# calculating bray curtis distance amongst samples from same age and profiled 
# with same method
amp_under15 <- calc_bc_dist(5, total_columns, "less than 15 months", "amp")
mgx_under15 <- calc_bc_dist(5, total_columns, "less than 15 months", "mgx")
amp_15to30 <- calc_bc_dist(5, total_columns, "15 to 30 months", "amp")
mgx_15to30 <- calc_bc_dist(5, total_columns, "15 to 30 months", "mgx")
amp_over30 <- calc_bc_dist(5, total_columns, "older than 30 months", "amp")
mgx_over30 <- calc_bc_dist(5, total_columns, "older than 30 months", "mgx")

bc_method_df <- rbind(amp_under15, mgx_under15, amp_15to30, mgx_15to30, 
                      amp_over30, mgx_over30)

bc_method_df$dev_stage <- factor(bc_method_df$dev_stage,
                                 levels = c("less than 15 months", 
                                            "15 to 30 months", 
                                            "older than 30 months"),
                                 ordered = TRUE)

# visualizing bray curtis distances
p2 <- ggplot(bc_method_df, aes(x = dev_stage, y = BC_dist)) + 
  geom_boxplot(aes(colour = method)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("Bray Curtis dissimilarity") + xlab("developmental stage") +
  labs(title = "", tag = "B") +
  theme(legend.position = c(0.9, 0.95)) + theme(legend.title=element_blank())
p2

# statistics

# does microbial community structure vary with profiling method?
my_dist <- vegdist(abund_table, "bray", diag=FALSE,upper=FALSE)
permanova <- adonis(my_dist ~ method, 
                      data=df, 
                      permutations=9999)


# paired vs. unpaired Bray-curtis 
bc_matrix <- as.matrix(vegdist(abund_table, "bray", diag=FALSE,upper=FALSE))

# upper left side of matrix is bc differences for mgx samples: mgx samples
mgx_bc <- bc_matrix[1:total_samples/2,1:total_samples/2]
mgx_bc[lower.tri(mgx_bc)] <- NA
mgx_bc[mgx_bc == 0.0] <- NA
mgx_bc_vec <- na.omit(as.vector(mgx_bc))
mgx_bc_vec

# lower right side of matrix = BC for 16S samples: 16S samples

amp_bc <- bc_matrix[((total_samples/2)+1):total_samples,((total_samples/2)+1):total_samples]
amp_bc[lower.tri(mgx_bc)] <- NA
amp_bc[amp_bc == 0.0] <- NA
amp_bc_vec <- na.omit(as.vector(amp_bc))
amp_bc_vec

# diagonal of the upper right side of matrix = BC for same samples (16S vs mgx)
paired_bc_vec <- as.vector(diag(bc_matrix[1:total_samples/2, 
                                          ((total_samples/2)+1):total_samples]))
length(paired_bc_vec) <- length(total_samples**2)

# samples with largest differences
paired_distances <- data.frame(distances = (na.omit(paired_bc_vec)), 
                               sampleid = as.vector(unique(df$sampleid)))
paired_distances$distances <- as.numeric(paired_distances$distances)

# upper right side of matrix excluding diagonal = BC for different samples 
# (16S vs mgx)
unpaired_bc <- diag.remove(bc_matrix[1:(total_samples/2), 
                                     ((total_samples/2)+1):total_samples], 
                           remove.val=NA)
unpaired_bc[lower.tri(unpaired_bc)] <- NA
unpaired_bc_vec <- na.omit(as.vector(unpaired_bc))

bc_df <- data.frame(paired_bc_vec, unpaired_bc_vec)

# for paper, add color to first two boxplots in figure
p0 <- ggplot(melt(bc_df), aes(variable, value)) + geom_boxplot()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("Bray Curtis dissimilarity") + xlab("between sample beta diversity")+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("16S & mgx, paired", "16S & mgx, unpaired"))+
  labs(title = "", tag = "C")
p0

# statistics
t.test(bc_df$paired_bc_vec, bc_df$unpaired_bc_vec)
mean(bc_df$unpaired_bc_vec) - mean(bc_df$paired_bc_vec, na.rm = TRUE)

mean(bc_df$paired_bc_vec, na.rm = TRUE)

# Bray-curtis for different ages based on profiling method
mean(bc_method_df[bc_method_df$dev_stage == "less than 15 months",]$BC_dist, 
     na.rm=TRUE)
mean(bc_method_df[bc_method_df$dev_stage == "15 to 30 months",]$BC_dist, 
     na.rm=TRUE)
mean(bc_method_df[bc_method_df$dev_stage == "older than 30 months",]$BC_dist, 
     na.rm=TRUE)


# calculating paired & unpaired Bray Curtis distances for kids of different ages
### see above function for arguments

calc_paired_BC <- function(bug_start, bug_end, dev_stage_arg){
  
  df_amp <- df[df$dev_stage == dev_stage_arg
               & df$method == "amp",]
  df_mgx <- df[df$dev_stage == dev_stage_arg
               & df$method == "mgx",]
  amp_abund <- df_amp[,bug_start:bug_end]
  mgx_abund <- df_mgx[,bug_start:bug_end]
  
  n_samples <- nrow(df_amp)
  
  combined_abund <- smartbind(amp_abund, mgx_abund)
  combined_abund[is.na(combined_abund)] <- 0
  combined_abund_matrix <- as.matrix(vegdist(combined_abund, "bray", 
                                             diag=FALSE,upper=FALSE))
  combined_abund_matrix[combined_abund_matrix == 0.0] <- NA
  
  paired <- as.vector(diag(combined_abund_matrix[1:n_samples, 
                                                 (n_samples+1):(2*n_samples)]))
  unpaired <- diag.remove(combined_abund_matrix[1:n_samples, 
                                                (n_samples+1):(2*n_samples)], 
                          remove.val=NA)
  unpaired[lower.tri(unpaired)] <- NA
  unpaired <- na.omit(as.vector(unpaired))
  
  paired_df <- data.frame(paired)
  colnames(paired_df)<-c("BC_dist")
  paired_df["method"] <- "paired"
  paired_df["dev_stage"] <- dev_stage_arg
  
  unpaired_df <- data.frame(unpaired)
  colnames(unpaired_df)<-c("BC_dist")
  unpaired_df["method"] <- "unpaired"
  unpaired_df["dev_stage"] <- dev_stage_arg
  return(smartbind(paired_df, unpaired_df))
}

paired_under15 <- calc_paired_BC(5, total_columns, "less than 15 months")
paired_15to30 <- calc_paired_BC(5, total_columns, "15 to 30 months")
paired_over30 <- calc_paired_BC(5, total_columns, "older than 30 months")

age_bc_df <- rbind(paired_under15, paired_15to30, paired_over30)
age_bc_df$dev_stage <- factor(age_bc_df$dev_stage,
                              levels = c("less than 15 months", 
                                         "15 to 30 months", 
                                         "older than 30 months"),ordered = TRUE)

# plotting paired bray curtis distances, broken down by age and method
p3 <- ggplot(age_bc_df, aes(dev_stage, BC_dist)) + geom_boxplot(aes(fill=method))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("Bray Curtis dissimilarity") + xlab("between sample beta diversity") +
  scale_fill_manual(values=c("#69b3a2", "grey")) + 
  theme(legend.title=element_blank()) + theme(legend.position = c(0.9, 0.95)) + 
  labs(title = "", tag = "C")
p3

# statistics
mean(age_bc_df[age_bc_df$dev_stage == "less than 15 months" & 
                 age_bc_df$method == "paired",]$BC_dist, 
     na.rm=TRUE)
mean(age_bc_df[age_bc_df$dev_stage == "15 to 30 months"& 
                 age_bc_df$method == "paired",]$BC_dist, 
     na.rm=TRUE)
mean(age_bc_df[age_bc_df$dev_stage == "older than 30 months" &
                 age_bc_df$method == "paired",]$BC_dist, 
     na.rm=TRUE)

t.test((age_bc_df[age_bc_df$dev_stage == "less than 15 months" & 
                    age_bc_df$method == "paired",]$BC_dist), 
       (age_bc_df[age_bc_df$dev_stage == "15 to 30 months" |
                    age_bc_df$dev_stage == "older than 30 months" & 
                    age_bc_df$method == "paired",]$BC_dist))

# making a RDA plot 

abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- df[,c("sampleid","AgeMonths","dev_stage","method")]

### Pcoa plot
dist <- vegdist(abund_table,  method = "bray")
PCOA <- pcoa(dist)

pcoa_vectors <- as.data.frame(PCOA$vectors)

pcoa_2_vectors <- cbind.data.frame(pcoa_vectors$Axis.1, pcoa_vectors$Axis.2)
pcoa_plot <- cbind(pcoa_2_vectors, meta_table$sampleid, meta_table$AgeMonths, 
                   meta_table$dev_stage, meta_table$method)

colnames(pcoa_plot)<-c("x","y","sampleid", "AgeMonths", "dev_stage", "method")
pcoa_plot$dev_stage <- factor(df_sites$dev_stage,
                              levels = c("less than 15 months", 
                                         "15 to 30 months", 
                                         "older than 30 months"),ordered = TRUE)

axis1 <- paste("PCoA 1, ", 
               round(PCOA$values$Relative_eig[1],4)*100, "%", sep = "")
axis2 <- paste("PCoA 2, ", 
               round(PCOA$values$Relative_eig[2],4)*100, "%", sep = "")

pcoa1 <- ggplot(data = pcoa_plot, aes(x,y,colour=dev_stage, shape = method, 
                                      group = sampleid))
pcoa1 <- pcoa1+geom_point(aes(colour=dev_stage, shape = method, 
                              group = sampleid), size = 3, alpha = 0.7) + 
  geom_line(size = 0.5) + theme_bw()+labs(x=axis1, 
                                          y=axis2, 
                                          color="developmental stage", 
                                          shape = "profiling method") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(legend.position = c(0.9, 0.99),
        legend.text=element_text(size=10),
        legend.title=element_blank()) + 
  labs(title = "", tag = "D")+
  theme(legend.key.size = unit(0.35, "cm"))
pcoa1


# variance explained by first 10 principal components
top_10_components <- (cumsum(PCOA$values$Relative_eig[1:10]))
names(top_10_components) <- seq(1, 10, by=1)


barplot(top_10_components,
        ylab="cummulative percent explained (%)",
        xlab="principal component #")


# aggregating all the plots together

gl <- list(p1, p2, p3, pcoa1)
grid.arrange(grobs = gl,widths = c(2, 1, 1), 
             layout_matrix = rbind(c(1,2, 2), c(3, 4, 4)))
