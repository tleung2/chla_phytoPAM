  ### Using chl_extractions dataset
  ### Using PhytoPAM_chla dataset

library(tidyverse)
library(Hmisc)
library(ggpubr)
library(ggpmisc)

#######################################################################
   ##################   IMPORT AND TIDY DATA   ###################

  ### Import data as chl_extractions as chla
chla=read.csv("chl_extractions.csv",header=T)

  ### remove sample B132 from chla dataframe
chla2<- subset(chla, sample_ID!= "B132")

  ### Import PhytoPAM_chla as PAMdata (using GUI to import excel)
  ### omit any NAs in PAM dataset
PAMdata2<- as.data.frame(na.omit(PAMdata))

PAMdata2$Month<-factor(PAMdata2$Month, levels = c("May", "June", "July",
                                                  "August", "September"))
  
  ### remove samples with fit error = 32000 (error in the field)
PAMdata3 <-subset(PAMdata2, Fit_error != "32000")

  ### Convert Week to factor
PAMdata3$Week <-factor(PAMdata3$Week)
   
  ### assign factor and levels for fit error
chla$fit_error <- factor(chla$fit_error, levels = c("0","0-1", "> 1"))     

#######################################################################
  #####################    SUBSETTING DATA    #####################

  ### Subset chla by PhytoPAM_chla range: 50, 100, 200
  ### Want to see if linearity is better at lower or higher PhytoPAM_chla
less50<-subset(chla2, PhytoPAM_chla <= 50)
less100<-subset(chla2, PhytoPAM_chla <= 100)
less200<-subset(chla2, PhytoPAM_chla <= 200)
greater200<-subset(chla2, PhytoPAM_chla >= 200)


#######################################################################
  ###################    SIGNIFICANCE TESTING   ##################
  
  ## Run Kruskal test
  ## if p < 0.05, then there is a significant difference
set.seed(123)
kruskal.test(fe.data3$Fe, fe.data3$Landform,correct=FALSE, na.rm = TRUE)

  ## Kruskal test show significant difference in Fe
  ## Run Wilcoxon pair test to see where difference is
pairwise.wilcox.test(fe.data3$Fe, fe.data3$site,
                     conf.int = 0.99, correct=FALSE, na.rm = TRUE)


#######################################################################
  ########################    BOXPLOTS   #########################

  ### --------   Boxplot: taxa specific chla each month   ----------

  ### Looks at range in chla for each taxa during each month
  ### 1) Pivot longer and re-shape data
  ### Make column for taxa and another for chlorophyll value
PAMdata4<- pivot_longer(PAMdata3, cols = c(8:11),
                        names_to = c("taxa"),
                        values_to = "chla")

  ### 2) Make boxplot and group by taxa
  ### scale_fill_discrete(labels = c()) --> manually assign legend labels
  ### make new column (order) so that plot is in order: blue, brown, green, PE 
PAMdata4$order <- factor(PAMdata4$taxa, levels = c("cyano_chla", "brown_chla",
                      "green_chla", "PE_chla"))
  ### rename months to make them shorter
PAMdata4$Month <- gsub("August", "Aug", PAMdata4$Month)
PAMdata4$Month <- gsub("September", "Sept", PAMdata4$Month)
PAMdata4$Month <- gsub("June", "Jun", PAMdata4$Month)
PAMdata4$Month <- gsub("July", "Jul", PAMdata4$Month)
PAMdata4$Month <- factor(PAMdata4$Month, levels = c("May", "Jun",
                                                   "Jul", "Aug", "Sept"))

  ### Make boxplot using ggplot
ggplot(PAMdata4, aes(x = Month, y = chla, fill = order)) +
  geom_boxplot() +
  labs(y = expression(paste('PhytoPAM Chl a (', mu, 'g/L)'))) +
  scale_fill_manual(values = c("#66CCCC", "#FFFF00", "#99CC00", "#FF9900"),
                    labels = c("'Blue' group", "'Brown' group",
                               "'Green' group","'Red' group")) +
  facet_wrap(.~order, scale = "free_y", ncol = 2) +
  theme(panel.background = element_blank(),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.line = element_line(color = "black"),
        strip.text = element_text(size = 14, color = "black"),
        legend.position = "none",
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(color = "white"))

  ### ------------  Boxplot: percent error  ----------------
  ### Shows if a specific lake consistently has high differences
  ### between phytoPAM and extracted chla
ggplot(chla, aes(x=Site, y= error)) + 
  geom_boxplot() +
  ylab("Percent Error") +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 24, color = "black"),
        panel.grid.major = element_line(color = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24, color = "black", 
                                   angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 24, color = "black"), ##text of graph titles in facet wrap
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        legend.position = "bottom",
        legend.key = element_blank())



#########################################################################
  ####################### LINEAR REGRESSION ##########################

  ### Correlate between PhytoPAM and extracted chla
  ### stat_smooth() shows the confidence interval
  ### stat_poly() calculates the linear equation
  ### stat_corr() calculated the correlation coefficient
ggplot(chla2, aes(x=Acetone_chla, y=PhytoPAM_chla)) + 
  geom_point(colour="black", size = 4) + 
  stat_smooth(method = 'lm', aes(color = 'linear'), se = TRUE) + ## Turns on confidence intervals
  stat_poly_eq(aes(label = ..eq.label..), formula = y ~ x, parse = TRUE, size = 6) +                                 ## Turns on equation
  stat_cor(label.x.npc = "center", label.y.npc = "top", size = 6) + ## Turns on r value
  labs(x = expression(paste('PhytoPAM-derived Chl a (', mu, 'g/L)')),
       y = expression(paste('Spectrophotometric Chl a (', mu, 'g/L)'))) +
  #scale_y_continuous(position = "right") +  ## places y scale on right
  facet_wrap(~fit_error, scales = "free", ncol = 3) +
  theme_classic() +
  theme(axis.text.y.left = element_text(size=28, color = "black"), 
        axis.text.x.bottom = element_text(size=28, color = "black"),
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28),
        strip.text = element_text(size = 28))

chla2<- chla %>%
  filter(Site %in% c("Prairie Rose", "Lake Three Fires", "Honey Creek Resort",
                     "Green Valley", "Lake Darling", "Union Grove"))

chla3<-chla %>%
  filter(Site %in% c("Lake Macbride", "Lake Keomah", "Big Creek",
                     "North Twin West", "North Twin East"))
#########################################################################
  ###################### STACKED AREA PLOT ###########################

    ## Transform phytoPAM chl-a columns for stacked area plotting
data3<-pivot_longer(pam_data,8:11,names_to = "phyto", values_to = "chl")

  ## Subset data by year
data.2018 <-subset(data3, year == "2018")
data.2018$phyto <- factor(data.2018$phyto,
                           levels=(c("PE-type_chla","brown_chla","green_chla","cyano_chla")))
viking<-subset(data.2020, Site == "Viking Lake")

ggplot(data.2020, aes(x = date, y = chl, fill = phyto), na.rm = TRUE) + 
  geom_area(aes(y = chl)) +
  #scale_fill_viridis(discrete = TRUE) +
  labs(y = expression(paste('Chlorophyll a (', mu, 'g/L)'))) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        panel.grid.major = element_line(color = "black"),
        axis.text = element_text(size = 18, color = "black"), 
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.line.y = element_line(color = "black")) + 
  facet_wrap(.~Site, scales = "free_y", ncol = 4) +
  scale_fill_manual(values=c("#CCCC00","#996600","#006600","#66CCCC"))
  

