
setwd("/Users/sophiafatima/Library/CloudStorage/OneDrive - MiddleburyCollege/Classes/THESIS/Hg Analysis")
load("Hg.RData")

hg_full<-read.csv(file="hg_raw.csv",  header=TRUE, sep=",")

library(dplyr)
library(ggplot2)
packageVersion("dplyr")
packageVersion("stats")
citation("stats")

hg_full<-hg_full%>% slice(1:123, 163:200)%>% # I pulled the samples in the run that had the very high concentrations and high cs456 at the end
  mutate(SampleName = case_match(
    SampleName, #combining values that are the same thing into one value
    c("acs_water","acswater", "acswater_real", "blank","proced_blank" ) ~ "proc_blank",
    c("water", "water blank")~ "water_blank", c("cs-456", "cs456")~"cs456", .default = SampleName))

#Collapse replicates into means for each replicate

hg<-hg_full%>% #plus removing autoblanks and spikes, and some other things
  filter(!grepl("auto", SampleName, ignore.case = TRUE))%>%
  filter(!grepl("spike", SampleName, ignore.case = TRUE))%>%
  filter(!SampleName %in% c( "procb", "ib", "clean"))%>%
  group_by(SampleName)%>%
  summarize(samplesize=n(), 
            mercury=mean(Concentration), 
            stdev=sd(Concentration))

#Limit of detection mean blank plus 3x stdev
# limit of quantification mean blank plus 10x stdev
hg%>%filter(SampleName=="water_blank")%>%summarize(LOD=mercury+stdev*3, LOQ=mercury+stdev*10)

#LOD   LOQ
#<dbl> <dbl>
#0.436  1.16

hg<-hg%>%mutate(LOQ_flag=ifelse(mercury>1.16, "NA", "flag"))
#all above the level of quantification


#your need to make one or two categorical variables to explain your sites
# and probably collapse data again by core site if you want to average all three cores

#download the csv of hg and add column w/river and site names
write.csv(hg, "C:/Users/hamia/OneDrive - Middlebury College/Classes/THESIS/Hg Analysis/hg.csv", row.names=TRUE)

hg_collapse<-read.csv(file="hg.csv",  header=TRUE, sep=",")

hg_collapse2<-hg_collapse%>%filter(River=="E" | River=="W")

write.csv(hg_collapse2, "/Users/sophiafatima/OneDrive - Middlebury College/Classes/THESIS/Hg Analysis/hg_collapse.csv", row.names=TRUE)

hg_collapse2%>% group_by(Type)%>%summarize(mean=mean(mercury), stdev=sd(mercury))

#Type   mean stdev
#<chr> <dbl> <dbl>
#1 B      9.52  3.42
#2 RL    25.7  25.8 
#3 RW    12.5  10.1 

#Levene's test for homogeneity of variance
library(car)
packageVersion("car")
leveneTest(mercury~River,data=hg_collapse2) #P=0.04073, so variance is NOT equal

#let's do a t-test to compare means between two rivers, E & W
#w/o taking into account the different “sites”

t.test(mercury~River,var.equal=F,data=hg_collapse2) #p-value = 0.004905

#rivers are significantly different, so I will analyze per river

### Englesby ###

#we first construct a linear model (lm) object
hg_E <- filter(hg_collapse2, River=="E")

hg_E%>% group_by(Type)%>%summarize(mean=mean(mercury), stdev=sd(mercury))

#Type   mean stdev
#<chr> <dbl> <dbl>
#  1 B      8.28  2.46
#2 RL    41.9  25.8 
#3 RW    16.3  12.2 

lm.E<-lm(mercury~Type,data=hg_E)

#next we extract the residuals of the lm object
#we will check that these are normally distributed
#residual takes each value of each sample and subtracts it from that group's mean
residuals<-residuals(lm.E)
shapiro.test(residuals)
#the p-value is 0.01952, so the data DO NOT follow a normal distribution

hist(residuals, data=hg_E)
#data is neither right-skewed nor left skewed I don't think taking log would work
#but let's try

lm.2<-lm(log(mercury)~Type, data=hg_E)
residuals2<-residuals(lm.2)
shapiro.test(residuals2) #p=0.883, so the data now DO follow a normal distribution

#we can also test for homogeneity of variance
leveneTest(mercury~Type,data=hg_E)
#p=0.0407, reject null hypothesis -> variances are NOT equal

#Let's do non-parametric analysis of variance

kruskal.test(mercury ~ Type,
             data = hg_E)

#	Kruskal-Wallis rank sum test

#data:  mercury by Type
#Kruskal-Wallis chi-squared = 12.436, df = 2, p-value = 0.001993

install.packages("FSA")
library(FSA)

dunnTest(mercury ~ Type,
         data = hg_E,
         method = "holm"
)

#Comparison         Z     P.unadj       P.adj
#1     B - RL -3.303014 0.000956516 0.002869548 **
#2     B - RW -1.395239 0.162943875 0.162943875
#3    RL - RW  2.431946 0.015017949 0.030035898 *

### Winooski ###
#we first construct a linear model (lm) object
hg_W <- filter(hg_collapse2, River=="W")

hg_W%>% group_by(Type)%>%summarize(mean=mean(mercury), stdev=sd(mercury))

#Type   mean stdev
#<chr> <dbl> <dbl>
#1 B     11.2   4.36
#2 RL     9.40 12.2 
#3 RW     8.71  5.75

lm.W<-lm(mercury~Type,data=hg_W)

#next we extract the residuals of the lm object
#we will check that these are normally distributed
#residual takes each value of each sample and subtracts it from that group's mean
residuals<-residuals(lm.W)
shapiro.test(residuals)
#the p-value is 0.0005769, so the data DO NOT follow a normal distribution

hist(residuals, data=hg_W)

#It looks like distribution is right-skewed so we take log
lm.2<-lm(log(mercury)~Site, data=cn_W)
residuals2<-residuals(lm.2)
shapiro.test(residuals2) #p=0.883, so the data now DO follow a normal distribution

#we can also test for homogeneity of variance
leveneTest(mercury~Type,data=hg_W)
#p=0.6829, accept null hypothesis -> variances ARE equal

#our assumptions are met, we can proceed to ANOVA
#the anova function requires an lm object
anova(lm.2)

#Analysis of Variance Table for Winooski

#Response: log(mercury)
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Type       2 7.3083  3.6541  11.457 0.0005442 ***
#  Residuals 19 6.0601  0.3190  

# TUKEY TEST
library(emmeans)
countpairs <- emmeans(lm.2, "Type")
countUnplanned <- contrast(countpairs, method = "pairwise", adjust = "tukey")
countUnplanned

#contrast estimate    SE df t.ratio p.value
#B - RL     -1.478 0.339 19  -4.356  0.0009 **
#B - RW     -0.537 0.339 19  -1.584  0.2767 
#RL - RW     0.941 0.266 19   3.534  0.0060 **

# you can then pull data from both the first collapse of replicates to make a geom_point, 
#and then another geom point pulling the average from the overall site


library(ggplot2)
  ggplot() + geom_point(data=hg_collapse, aes(x = Type, y = mercury), show.legend=FALSE)+
  labs(y = "Total Hg ng/g", x=NULL) + 
  theme_bw() + theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust=1 ))
  
#I took out below line of code from second line of above chunk:
#  geom_point(data=sophia_collapse, aes(x=location, y=mean))+
  
  ggplot() + geom_point(data=hg_E, aes(x = Type, y = mercury, color=Site), show.legend=FALSE)+
  labs(y = "Total Hg ng/g", x=NULL) + 
  theme_bw() + theme(panel.border = element_blank())+
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust=1 ))
  
  ggplot() + geom_point(data=hg_W, aes(x = Type, y = mercury, color=Site), show.legend=FALSE)+
    labs(y = "Total Hg ng/g", x=NULL) + 
    theme_bw() + theme(panel.border = element_blank())+
    theme(text=element_text(size=20),
          axis.text.x = element_text(angle = 45, hjust=1 ))
  
  #making a box-plot
  hg_E%>%
    ggplot(aes(x=Type, y=mercury, color))+
    geom_boxplot()+
    xlab("Englesby Sites")+
    ylab("Mercury (ng/g)")+
    theme_classic()
  
  hg_W%>%
    ggplot(aes(x=Type, y=mercury))+
    geom_boxplot()+
    xlab("Winooski Sites")+
    ylab("Mercury (ng/g)")+
    theme_classic()
 

  write.csv(hg_collapse, "hg_data.csv")

  
#### REDO WITH AVERAGED MERCURY VALUES ####
  
hg_avg<-hg_collapse2%>%
    group_by(Site)%>%
    mutate(Average_Hg=mean(mercury),
          stdev_Average_Hg=sd(mercury)
          )%>%
    select(Site, River, Type, Average_Hg, stdev_Average_Hg)%>%
    distinct()
  
#Levene's test for homogeneity of variance
  library(car)
leveneTest(Average_Hg~River,data=hg_avg) #P=0.1364, so variance IS equal
  
  #let's do a t-test to compare means between two rivers, E & W
  #w/o taking into account the different “sites”
  
  t.test(Average_Hg~River,var.equal=F,data=hg_avg) #p-value = 0.07116
  #No significant difference in Hg between rivers??
  
#Welch Two Sample t-test
  
#  data:  Average_Hg by River
#  t = 2.0975, df = 7.5711, p-value = 0.07116
#  alternative hypothesis: true difference in means between group E and group W is not equal to 0
#  95 percent confidence interval:
#    -1.853425 35.443964
#  sample estimates:
#    mean in group E mean in group W 
#  26.150095        9.354825 
  
#CONTINUE WITH Hg ANALYSIS BY SITE FOR ALL SAMPLES
  
hg_avg%>% group_by(Type)%>%summarize(mean=mean(Average_Hg), stdev=sd(stdev_Average_Hg))

# Type   mean stdev
# <chr> <dbl> <dbl>
#1 B      9.72  1.34
#2 RL    25.7  12.6 
#3 RW    12.5   6.62  

lm.all<-lm(Average_Hg~Type,data=hg_avg)

residuals<-residuals(lm.all)
shapiro.test(residuals)
#the p-value is 0.1587, so the data DO follow a normal distribution

#we can also test for homogeneity of variance
leveneTest(Average_Hg~Type,data=hg_avg)
#p=0.1086, reject null hypothesis -> variances ARE equal

anova(lm.all)

#Analysis of Variance Table

#Response: Average_Hg
#           Df  Sum Sq Mean Sq F value Pr(>F)
#Type       2  668.79  334.39  1.2215 0.3319
#Residuals 11 3011.45  273.77   

###CONCLUSION###
#No significant difference in mercury level across sites

###SUBSETTING E###
    
hg_E_avg <- filter(hg_avg, River=="E")

hg_E_avg%>% group_by(Type)%>%summarize(mean=mean(Average_Hg), stdev=sd(stdev_Average_Hg))

#Type   mean stdev
#<chr> <dbl> <dbl>
#  1 B      8.28 NA   
#2 RL    41.9  15.8 
#3 RW    16.3   8.53

lm.E<-lm(Average_Hg~Type,data=hg_E_avg)

#next we extract the residuals of the lm object
#we will check that these are normally distributed
#residual takes each value of each sample and subtracts it from that group's mean
residuals<-residuals(lm.E)
shapiro.test(residuals)
#the p-value is 0.2312, so the data DO follow a normal distribution

#we can also test for homogeneity of variance
leveneTest(Average_Hg~Type,data=hg_E_avg)
#p=0.7774, reject null hypothesis -> variances ARE equal

#our assumptions are met, we can proceed to ANOVA
#the anova function requires an lm object
anova(lm.E)

#Analysis of Variance Table

#Response: Average_Hg
#           Df Sum Sq Mean Sq F value Pr(>F)
#Type       2 1355.9  677.95  2.6576 0.1844
#Residuals  4 1020.4  255.10 

### SUBSETTING W

hg_W_avg <- filter(hg_avg, River=="W")

hg_W_avg%>% group_by(Type)%>%summarize(mean=mean(Average_Hg), stdev=sd(stdev_Average_Hg))

#Type   mean  stdev
#<chr> <dbl>  <dbl>
#  1 B     11.2  NA    
#2 RL     9.40  8.00 
#3 RW     8.71  0.507

lm.W<-lm(Average_Hg~Type,data=hg_W_avg)

#next we extract the residuals of the lm object
#we will check that these are normally distributed
#residual takes each value of each sample and subtracts it from that group's mean
residuals<-residuals(lm.W)
shapiro.test(residuals)
#the p-value is 0.139, so the data DO follow a normal distribution

#we can also test for homogeneity of variance
leveneTest(Average_Hg~Type,data=hg_W_avg)
#p=0.807, reject null hypothesis -> variances ARE equal

anova(lm.W)

#Analysis of Variance Table

#Response: Average_Hg
#           Df  Sum Sq Mean Sq F value Pr(>F)
#Type       2   4.492   2.246  0.0288 0.9718
#Residuals  4 312.151  78.038   

write.csv(hg_avg, "/Users/sophiafatima/Library/CloudStorage/OneDrive-MiddleburyCollege/Classes/THESIS/Hg Analysis/hg_avg.csv", row.names=TRUE)

#Take average of each type of sample: WRW, WRL, WB, ERL, ERW, EB
hg_avg_final<-hg_avg%>%
  group_by(River, Type)%>%
  mutate(Average_Hg_samplingType=mean(Average_Hg),
         se_Average_Hg_samplingType=sd(Average_Hg)/sqrt(length(Average_Hg)) #get standard error
  )%>%
  # Replace NA standard errors using stdev_Average_Hg and hardcoded count
  ungroup() %>%
  mutate(
    se_Average_Hg_samplingType = case_when(
      Site == "eb" & is.na(se_Average_Hg_samplingType) ~ stdev_Average_Hg / sqrt(4),
      Site == "wb" & is.na(se_Average_Hg_samplingType) ~ stdev_Average_Hg / sqrt(3),
      TRUE ~ se_Average_Hg_samplingType
    )
  ) %>%
  select(River, Type, Average_Hg_samplingType, Average_Hg, se_Average_Hg_samplingType) %>%
  distinct()

hg_avg_final<-hg_avg_final%>%
  mutate(sampling_site_type = paste(River, Type, sep = ""))

###SCATTERPLOT

#Making scatterplot with Winooski data
hg.plot<-ggplot() + geom_point(data=hg_avg_final, aes(x = sampling_site_type, y = Average_Hg, color=River, shape=Type), size=5, show.legend=TRUE)+
  # Group averages (from Average_Hg_samplingType column)
  geom_point(data = hg_avg_final, 
             aes(x = sampling_site_type, y = Average_Hg_samplingType, color = River), 
             shape = 18, size = 5, stroke = 1.5, show.legend = FALSE) +
  scale_shape_manual(values=c(0, 1, 2, 4, 5, 6, 7))+
  labs(y = "Average THg(ng/g)", x=NULL) + 
  theme_bw() + theme(panel.border = element_blank())+
  theme(text=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust=1 ))+
  labs(color="Stream", shape="Site", size=NULL)

#Prettification of the plot above:

hg.plot <- ggplot() + 
  
  # Individual sample points
  geom_point(data = hg_avg_final, 
             aes(x = sampling_site_type, y = Average_Hg, color = River, shape = Type), 
             size = 4, alpha = 1, stroke = 0.7) +
  
  # Group means (Average_Hg_samplingType)
  geom_point(data = hg_avg_final, 
             aes(x = sampling_site_type, y = Average_Hg_samplingType, color = River),
             shape = 18, size = 6, show.legend = FALSE) +
  
  # Axis & legend styling
  scale_shape_manual(values = c(0, 1, 2, 4, 5, 6, 7)) +
  
  labs(
    y = expression(paste("Total Hg (ng/g dry weight)")),
    x = NULL,
    color = "Stream",
    shape = "Site"
  ) +
  
  # Clean, classic theme
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave("hg.tiff", plot = hg.plot, dpi = 1000, width = 8, height = 5, units = "in")
ggsave("hg.png", plot = hg.plot, dpi = 1000, width = 8, height = 5, units = "in")

