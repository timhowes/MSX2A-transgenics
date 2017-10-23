## Experiment 1 - PAXB eggs injected with 5.6 kb MSX2A BAC fragment and lens GFP marker

library(ggplot2)
d <- read.csv("transient_experiment_1.csv", header=TRUE)


## Embryos were initially sorted based on whether or not fluorescence was observed in at least one eye.
## Different groups (uninjected, injected but not fluorescent, and transgenic/fluorescent) were raised in different tanks.
## GFP fluorescence in the lens was checked again in adults at time of measurement.
## Measurement were made for dorsal spine 1, dorsal spine 2, dorsal spine 3, anal spine,
##   lower jaw, pectoral fin (measured along the middle), caudal fin (measured along the middle),
##   left and right pelvic spines
## All lengths in millimeters

## dorsal spine 1 and pelvic spines were frequently absent, as is typical for the PAXB line
## dorsal spine 3 and anal spine were very small, not scored if too small for caliper measurement


## abbreviate some names in the table
names(d)[4] <- "gfp.adult"   # green.eye.or.eyes.observed.as.adult == "gfp.adult"
levels(d$group) <-  c("no.gfp", "gfp", "uninj")
# "no.gfp" ==  injected, no GFP observed in embryo
# "gfp"    ==  injected, GFP observed in embryo
# "uninj"  ==  uninjected sibling


## Compare length of dorsal spine 2 between gfp and uninjected groups

## Take residuals to standard length based on linear regression
lm.ds2 <- lm(d.spine.2 ~ std.length, d)
plot( d$std.length, d$d.spine.2)
abline(lm(d.spine.2 ~ std.length, data=d), lty=2)
rsd.ds2 <- d$d.spine.2 - predict(lm.ds2, d)      # store residuals in this vector
rsd.ds2.1 <- rsd.ds2  # store this version for later recall


## Compare gfp group to all others (no.gfp and uninj)
d$gfp.vs.other <- d$group == "gfp"
t.test( rsd.ds2 ~ gfp.vs.other, d)       #p=.03319   t(12.2) = 2.4004
wilcox.test( rsd.ds2 ~ gfp.vs.other, d)  #p=.04053   W=29
                                         #difference in mean = 0.20672967 - -0.08039487 = 0.2871245
                                         #number in each group:  gfp 7, other 18


## Spine length residuals: fish 2, 4, and 9 are outliers in control group
ggplot(d, aes(gfp.vs.other, rsd.ds2)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=d$fish) +
    scale_x_discrete(labels=c("Other (no gfp and uninj)", "GFP"))


## Spine length residuals by tank
d$tank <- as.factor(d$tank)
ggplot(d, aes(tank, rsd.ds2)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=d$fish)

## Standard lengths by tank
ggplot(d, aes(tank, std.length)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=d$fish)


## Try regression based only on uninjected fish to calculate residuals
qplot(std.length, d.spine.2, data=d, color=group)  # some gfp fish above diagonal, some uninjected below diagonal
plot(d$std.length, d$d.spine.2)
abline(lm.ds2)              # original regression line
lm.ds2.uninj <- lm(d.spine.2 ~ std.length, d[d$group=="uninj",])
abline(lm.ds2.uninj, lty=2) # fitting only to uninjected (dashed line)
                            # seems less appropriate, as regression line
                            # moves away from the bulk of the data
rsd.ds2 <- d$d.spine.2 - predict( lm.ds2.uninj, d )
t.test( rsd.ds2 ~ gfp.vs.other, d )      #p=.04623   t(11.261) = 2.2392
wilcox.test( rsd.ds2 ~ gfp.vs.other, d)  #p=.07371   W=33
                                         #difference in mean = 0.4215429 - 0.1393920 = 0.2821509
                                         #number in each group:  gfp 7, other 18


## Try removing outliers from the control group
quantile(rsd.ds2.1[ !d$gfp.vs.other])
iqr <- IQR( rsd.ds2.1[ !d$gfp.vs.other])
q1 <- quantile( rsd.ds2.1[ !d$gfp.vs.other] )[2]
q3 <- quantile( rsd.ds2.1[ !d$gfp.vs.other] )[4]
r1 <- q1 - 1.5*iqr  #-.351707   residuals below this are outliers
r2 <- q3 + 1.5*iqr  #0.3029504  residuals above this are outliers

## Remove outliers, store data frame as d2
d2 <- d[-c(2,4,9),]

lm.ds2 <- lm(d.spine.2 ~ std.length, d2)
plot( d2$std.length, d2$d.spine.2)
abline(lm(d.spine.2 ~ std.length, data=d2), lty=2)
rsd.ds2.1b <- d2$d.spine.2 - predict(lm.ds2, d2)

t.test( rsd.ds2.1b ~ gfp.vs.other, d2)       #p=.04997   t(7.6853) = 2.3229
wilcox.test( rsd.ds2.1b ~ gfp.vs.other, d2)  #p=.04652   W=24
                                             #difference in mean = 0.15580820 - -0.07271049 = 0.2285187
                                             #number in each group:  gfp 7, other 15

ggplot(d2, aes(gfp.vs.other, rsd.ds2.1b)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=d2$fish) +
    scale_x_discrete(labels=c("Other (no gfp and uninj)", "GFP"))



library(pwr)
# Pooled standard deviation
sqrt((sd(rsd.ds2.1[d$gfp.vs.other])**2 + sd(rsd.ds2.1[!d$gfp.vs.other])**2)/2)   #.275

# Power calculation, t-test with unequal samples size
pwr.t2n.test(n1 = 7, n2 = 18, d = .5, sig.level=.05)   # power = 0.190, difference .14 mm
pwr.t2n.test(n1 = 7, n2 = 18, d = .8, sig.level=.05)   # power = 0.406, difference .22 mm
pwr.t2n.test(n1 = 7, n2 = 18, d = 1, sig.level=.05)    # power = 0.576, difference .28 mm
pwr.t2n.test(n1 = 7, n2 = 18, d = 1.3, sig.level=.05)  # power = 0.798, difference .36 mm








## Experiment 2 - PAXB eggs injected with 5.6 kb MSX2A BAC fragment and lens GFP marker

library(ggplot2)
m <- read.csv("transient_experiment_2.csv", header=TRUE)

## Remove fish 15. It exhibited tail deformities, only a small bit of GFP fluorescence
##   in left eye, and it had no dorsal spines at all
m <- m[-15,]



## abbreviate some names in the table
names(m)[5] <- "gfp.adult"   # green.eye.or.eyes.observed.as.adult == "gfp.adult"
levels(m$group) <-  c("no.gfp", "gfp", "uninj")
# "no.gfp" ==  injected, no GFP observed in embryo
# "gfp"    ==  injected, GFP observed in embryo
# "uninj"  ==  uninjected sibling


## Compare length of dorsal spine 2 between gfp and uninjected groups

## Take residuals to standard length based on linear regression
lm.ds2 <- lm(d.spine.2 ~ std.length, m)
rsd.ds2 <- m$d.spine.2 - predict(lm.ds2, m)
rsd.ds2.2 <- rsd.ds2  # store this version for later recall

## Compare gfp group to all others (no.gfp and uninj)
m$gfp.vs.other <- m$group == "gfp"
t.test( rsd.ds2 ~ gfp.vs.other, m)       #p=.08048   t(21.114) = 1.836
wilcox.test( rsd.ds2 ~ gfp.vs.other, m)  #p=.1104    W=230
                                         #difference in mean = 0.10886734 - -0.03313354 = 0.1420009
                                         #number in each group:  gfp 14, other 46
 
ggplot(m, aes(gfp.vs.other, rsd.ds2)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=m$fish) +
    scale_x_discrete(labels=c("Other (no gfp and uninj)", "GFP"))


## DS2 length residuals by tank
m$tank <- as.factor(m$tank)
ggplot(m, aes(tank, rsd.ds2)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=m$fish)

## Standard length by tank
## Uninjected tanks (1 and 5) are larger overall because fewer fish survived (7 per tank)
## Fish in GFP tank (2, 14 fish) are smaller on average, and those in
##  "injected, no gfp" tanks (3 and 4) are smallest (16 fish per tank)
ggplot(m, aes(tank, std.length)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=m$fish)

## T-test for standard length in uninjected tanks vs. others
m$uninjected.vs.other <- m$group == "uninj"
t.test( std.length ~ uninjected.vs.other, m )       #p=.0004117   t(56.62) = 3.7545
wilcox.test( std.length ~ uninjected.vs.other, m )  #p=.03571     W=202
                                                    #difference in mean = 64.38286 - 60.68783 = 3.69503
                                                    #number in each group:  uninj 14, other 46


## Compare gfp group to no.gfp injected group (ignore uninjected because of size differences)
m2 <- m[ m$group != "uninj", ]
lm.ds2.2b <- lm(d.spine.2 ~ std.length, m2)
rsd.ds2.2b <- m2$d.spine.2 - predict(lm.ds2.2b, m2)
t.test( rsd.ds2.2b ~ group, m2)       #p=.04906   t(25)=2.0687
wilcox.test( rsd.ds2.2b ~ group, m2)  #p=.07512   W=149
                                      #difference in mean = 0.11784210 - -0.05155592 = 0.169398
                                      #number in each group:  no.gfp 32, gfp 14

ggplot(m2, aes(group, rsd.ds2.2b)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=m2$fish) +
    scale_x_discrete(labels=c("Injected, no GFP","GFP"))



library(pwr)
# Pooled standard deviation
sqrt((sd(rsd.ds2.2[m$gfp.vs.other])**2 + sd(rsd.ds2.2[!m$gfp.vs.other])**2)/2)   #.252

# Power calculation, t-test with unequal samples size
pwr.t2n.test(n1 = 14, n2 = 46, d = .5, sig.level=.05)   # power = 0.364, difference .13 m
pwr.t2n.test(n1 = 14, n2 = 46, d = .8, sig.level=.05)   # power = 0.732, difference .20 mm
pwr.t2n.test(n1 = 14, n2 = 46, d = .9, sig.level=.05)   # power = 0.826, difference .23 mm
pwr.t2n.test(n1 = 14, n2 = 46, d = 1, sig.level=.05)    # power = 0.896, difference .25 mm
pwr.t2n.test(n1 = 14, n2 = 46, d = 1.3, sig.level=.05)  # power = 0.987, difference .33 mm








## Experiment 3 - PAXB cross from a stable line carrying 5xMSX2ACNE-eGFP-2A-MSX2A construct

## Fish carrying this construct have low-level GFP expression in eyes and in mouth, nostrils, pectoral fins
##   (spine GFP expression and median fin expression visible in early larval stages but difficult to see in adults)


library(ggplot2)
n <- read.csv("stable_line_experiment.csv")



## Compare length of dorsal spine 2 between transgenic and wild-type groups ("gfp" and "wt")

## Take residuals to standard length based on linear regression
lm.ds2 <- lm(d.spine.2 ~ std.length, n)
rsd.ds2 <- resid(lm.ds2)
rsd.ds2.3 <- rsd.ds2  # store this version for later recall


## Compare gfp group to wt siblings
t.test( rsd.ds2 ~ group, n)                  #p=.03403   t(22.602) = 2.2564
wilcox.test( rsd.ds2 ~ group=="gfp", n)      #p=.02916   W=55
                                      #difference in mean = 0.08043809 - -0.07507555 = 0.1555136
                                      #number in each group:  gfp 14, wt 15

ggplot(n, aes(group=="gfp", rsd.ds2)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=n$fish) +
    scale_x_discrete(labels=c("WT", "GFP"))



## Try regression based only on wt fish to calculate residuals
qplot(std.length, d.spine.2, data=n, color=group)
plot(n$std.length, n$d.spine.2)
abline(lm.ds2)
lm.ds2.wt <- lm(d.spine.2 ~ std.length, n[n$group=="wt",])
abline(lm.ds2.wt, lty=2)
rsd.ds2 <- n$d.spine.2 - predict( lm.ds2.wt, n )
t.test( rsd.ds2 ~ group, n)                       #p=.02713   t(22.198) = 2.3659
wilcox.test( rsd.ds2 ~ group=="gfp", n)           #p=.03277   W=56
                                        #difference in mean = 1.630145e-01 - 8.141654e-16 = 0.1630145
                                        #number in each group:  gfp 14, wt 15

ggplot(n, aes(group=="gfp", rsd.ds2)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=n$fish) +
    scale_x_discrete(labels=c("WT", "GFP"))





## Standard length by tank/group
t.test( std.length ~ tank, n)       # some smaller fish in wt tank
ggplot(n, aes(as.factor(tank), std.length)) + geom_boxplot(outlier.color="green") + geom_jitter(position=position_jitter(width=.14)) + geom_text(label=n$fish)

## Histogram comparison of gfp group and wt group
ggplot(n, aes(x=rsd.ds2,fill=group)) + geom_histogram(binwidth=.1, alpha=.5, position="identity")



## Compare gfp group to wt siblings, removing outliers
n2 <- n[ -c(6,13), ]
rsd.ds2.3b <- resid( lm(d.spine.2 ~ std.length, n2 ) )
t.test( rsd.ds2.3b ~ group, n2)                 #p=.00513   t(16.357) = 3.2289
wilcox.test( rsd.ds2.3b ~ group=="gfp", n2)     #p=.005596  W=35
                                      #difference in mean = 0.09514592 - -0.10246483 = 0.1976108
                                      #number in each group:  gfp 14, wt 15

ggplot(n2, aes(group=="gfp", rsd.ds2.3b)) + geom_boxplot(outlier.color="green") +
    geom_jitter(position=position_jitter(width=.14)) + geom_text(label=n2$fish) +
    scale_x_discrete(labels=c("WT", "GFP"))
# Fish 5 is now considered an outlier too, but close to threshold



library(pwr)
# Pooled standard deviation
sqrt((sd(rsd.ds2.3[n$group=="gfp"])**2 + sd(rsd.ds2.3[n$group=="wt"])**2)/2)   #.184

# Power calculation, t-test with unequal samples size
pwr.t2n.test(n1 = 14, n2 = 15, d = .5, sig.level=.05)   # power = 0.255, difference .09 mm
pwr.t2n.test(n1 = 14, n2 = 15, d = .8, sig.level=.05)   # power = 0.546, difference .15 mm
pwr.t2n.test(n1 = 14, n2 = 15, d = 1, sig.level=.05)    # power = 0.737, difference .18 mm
pwr.t2n.test(n1 = 14, n2 = 15, d = 1.1, sig.level=.05)  # power = 0.814, difference .20 mm
pwr.t2n.test(n1 = 14, n2 = 15, d = 1.3, sig.level=.05)  # power = 0.921, difference .24 mm

# Power vs. sample size at effect size d = .8
p.out <- pwr.t.test( d = .8, sig.level=.05, power=.8)
print(p.out)
plot(p.out)  # sample size 26

# Power vs. sample size at effect size d = 1
p.out <- pwr.t.test( d = 1, sig.level=.05, power=.8)
print(p.out)
plot(p.out)  # sample size 17










##Final plots for Experiments 1-3

## Transient transgenic injection 1 (7 transgenics, 18 controls)
pdf("transgenic_plot_1.pdf", width=4.044, height=4.8, family="Helvetica")
ggplot(d, aes(gfp.vs.other,rsd.ds2.1)) +
    geom_boxplot(aes(fill=gfp.vs.other),fill=NA,outlier.shape=NA,width=.3,show.legend=FALSE, color="black") +
    geom_jitter(aes(x=as.numeric(gfp.vs.other)+1.3),position=position_jitter(width = .04),size=2.7,shape=21) +
    theme(panel.background = element_rect(fill="white", color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y=element_text(color="black",size=14), axis.title.y=element_text(size=18),
          axis.text.x=element_text(color="black",size=18), axis.ticks=element_line(color="black"),
          plot.margin = unit(c(5,5,2.5,4), "mm"), axis.ticks.length=unit(4, "pt"),
          plot.background = element_rect(fill="transparent")) +
    scale_y_continuous(limit=c(-.9,.9), breaks=c(-.9,-.6,-.3,0,.3,.6,.9)) +
    scale_x_discrete(labels=c("Control","Transgenic")) + ylab("Dorsal spine 2 length residuals (mm)") + xlab(NULL)
dev.off()

## Transient transgenic injection 2 (14 transgenics, 46 controls)
pdf("transgenic_plot_2.pdf", width=4, height=4.8, family="Helvetica")
ggplot(m, aes(gfp.vs.other,rsd.ds2.2)) +
    geom_boxplot(aes(fill=gfp.vs.other),fill=NA,outlier.shape=NA,width=.3,show.legend=FALSE, color="black") +
    geom_jitter(aes(x=as.numeric(gfp.vs.other)+1.3),position=position_jitter(width = .04),size=2.7,shape=21) +
    theme(panel.background = element_rect(fill="white", color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y=element_text(color="black",size=14), axis.title.y=element_text(size=18),
          axis.text.x=element_text(color="black",size=18), axis.ticks=element_line(color="black"),
          plot.margin = unit(c(5,5,2.5,4), "mm"), axis.ticks.length=unit(4, "pt"),
          plot.background = element_rect(fill="transparent")) +
    scale_y_continuous(limit=c(-.603,.6)) +
    scale_x_discrete(labels=c("Control","Transgenic")) + ylab("   ") + xlab(NULL)
dev.off()

## Stable transgenic construct (14 transgenics, 15 controls)
pdf("transgenic_plot_3.pdf", width=4, height=4.8, family="Helvetica")
ggplot(n, aes(group=="gfp",rsd.ds2.3)) +
    geom_boxplot(aes(fill=group=="gfp"),fill=NA,outlier.shape=NA,width=.3,show.legend=FALSE, color="black") +
    geom_jitter(aes(x=(as.numeric(group=="gfp"))+1.3),position=position_jitter(width = .04),size=2.7,shape=21) +
    theme(panel.background = element_rect(fill="white", color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y=element_text(color="black",size=14), axis.title.y=element_text(size=18),
          axis.text.x=element_text(color="black",size=18), axis.ticks=element_line(color="black"),
          plot.margin = unit(c(5,5,2.5,4), "mm"), axis.ticks.length=unit(4, "pt"),
          plot.background = element_rect(fill="transparent")) +
    scale_y_continuous(limit=c(-.6,.6)) +
    scale_x_discrete(labels=c("Control","Transgenic")) + ylab("  ") + xlab(NULL)
dev.off()
