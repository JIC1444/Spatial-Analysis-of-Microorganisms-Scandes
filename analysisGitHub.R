# Overview ----------------------------------------------------------------

#Title: Driving Factors Behind Change In Bacterial, Fungal and Plant Communities
#Over An Elevational Gradient

# This r script sets out to load in the *private* data provided to me from an 
# ongoing study in the Scandes Mountains in Norway.
# I focus on the interpretability of the visualisations of said data, by use of 
# 3 figures with appropriate statistical tests: a Knowledge graph, a box plot 
# with an ANOVA test (as well as a TukeyHSD test) and finally a scatter plot of 
# a linear regression model, including  p values.

# Packages ----------------------------------------------------------------

install.packages(c("tidyr", "readr", "ggplot2", "glue", "dplyr", "devtools", 
                   "vegan", "reticulate"))

library(tidyr)
citation("tidyr")
#Wickham H, Vaughan D, Girlich M (2024). _tidyr: Tidy Messy Data_. R package
#version 1.3.1, <https://CRAN.R-project.org/package=tidyr>.

library(readr)
citation("readr")
#Wickham H, Hester J, Bryan J (2024). _readr: Read Rectangular Text Data_. R
#package version 2.1.5, <https://CRAN.R-project.org/package=readr>.

library(ggplot2)
citation("ggplot2")
#H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New
#York, 2016.

library(glue)
citation("glue")
#Hester J, Bryan J (2024). _glue: Interpreted String Literals_. R package
#version 1.7.0, <https://CRAN.R-project.org/package=glue>.

library(dplyr)
citation("dplyr")
# Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar
#of Data Manipulation_. R package version 1.1.4,
#<https://CRAN.R-project.org/package=dplyr>.

library(devtools)
citation("devtools")
#Wickham H, Hester J, Chang W, Bryan J (2022). _devtools: Tools to Make
#Developing R Packages Easier_. R package version 2.4.5,
#<https://CRAN.R-project.org/package=devtools>.

library(vegan)
citation("vegan")
#Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R,
#Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M, Bedward M, Bolker B,
#Borcard D, Carvalho G, Chirico M, De Caceres M, Durand S, Evangelista H,
#FitzJohn R, Friendly M, Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D,
#Ouellette M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J (2022).
#_vegan: Community Ecology Package_. R package version 2.6-4,
#<https://CRAN.R-project.org/package=vegan>.

library(reticulate)
citation("reticulate")
#Ushey K, Allaire J, Tang Y (2024). _reticulate: Interface to 'Python'_. R
#package version 1.35.0, <https://CRAN.R-project.org/package=reticulate>.


# I created the Py2RNXGraph package as a Python extension to use the library 
# NetworkX. IGraph for RStudio does not have the rich functionality that 
# NetworkX does, which is what drove my descision to do this

# Requires version 3.11.7 of Python, a Virtual environment as well as libraries: 
# Pandas, NetworkX, Plotly and Matplotlib. This is all done in the lines below.
install_github("JIC1444/py2RNetworkX")
library(Py2RNetworkX)

version <- "3.11.7" 
install_python(version=version)
virtualenv_create("python-session", python_version=version)
use_virtualenv("python-session",require=TRUE)

libraries <- list("pandas","networkx","plotly","matplotlib")
for (lib in libraries){
  virtualenv_install(envname="python-session", 
                     lib,
                     ignore_installed=FALSE, 
                     pip_options=character())
}


# Functions ---------------------------------------------------------------

loadDataToTxt<-function(fileName, fileLocation){
  fileName<-read_table(file=fileLocation, col_names=TRUE)
}


# Function groups rows in data frame together, with the use of the taxonomy file
# which contains information about the phylum, class, order, etc. of the microbes
groupBy<-function(countTable, taxonomicTable, taxonomicGrouping){
  groups<-unique(taxonomicTable[[taxonomicGrouping]])
  groups<-groups[-which(is.na(groups))]
  out<-sapply(groups, function(group){
    ids<-taxonomicTable$ASV[which(taxonomicTable[[taxonomicGrouping]]==group)]
    idNums=which(countTable[[1]] %in% ids)
    colSums(countTable[idNums,-1], na.rm=T)
  })
  out=out[,which(colSums(out, na.rm=T)>0)]
  out
}


# Function which draws a box plot, then does a two-way anova, then a TukeyHSD
# Boxplot used as they are good for comparing distributions of sub-groups 
# within data (eg. the sets of temperatures and precipitations in the data)

# Two way ANOVA test used to see whether there were interactions between 
# temperature and precipitation factors

# TukeyHSD used to test whether there is any statistically significant pairs in
# the data across the sub-groups
boxPlotTwoWayAnova<-function(phylumDF){
  
  diversityVal<-diversity(phylumDF, index="shannon")
  ids<-rownames(phylumDF)
  idNums<-which(about$site %in% ids)
  
  df<-data.frame(diversity=diversityVal, precip<-about[idNums,"precipitation"],
                 temp=about[idNums, "temperature"])
  g<-ggplot(df, 
            aes(x=as.factor(precipitation), 
                y=diversity,
                fill=as.factor(temperature)))+
            geom_boxplot()+ 
            labs(x = "Expected Summer Precipitation (mm)",
                 y = "Shannon Diversity Index",)+ 
            scale_fill_manual(values=c("tomato", "skyblue", "#F0E442"))+
            guides(fill=guide_legend(title="Expected\nSummer\nTemperature\n(C)"))+
            theme_classic(base_size=20)+
            expand_limits(x=0, y=0)
  print(g)
  
  res.aov<-aov(diversityVal~as.factor(precipitation)*as.factor(temperature),
                 data=df)
  print(summary(res.aov))
  
  res.lm<-lm(diversityVal~as.factor(precipitation)*as.factor(temperature),
               data = df)
  print(summary(res.lm))
  
  print(TukeyHSD(res.aov))
}


# Function outputs the p values from the linear regression with precipitation and 
# temperature as the input variables and abundance as the output variable
pvalsViaRegression<-function(df){
  pvals=apply(df, 2, function(class){
          
    df=data.frame(abundance=class, 
                  precip=about[idNums, "P.grid"], 
                  temp=about[idNums, "T.grid"])
                  
    out=lm(abundance ~ P.grid + T.grid, data=df)
    sumres=summary(out)
    pf(sumres$fstatistic[1L], 
       sumres$fstatistic[2L], 
       sumres$fstatistic[3L], 
       lower.tail=FALSE)
  })
  print(p.adjust(pvals))
}


# Function plots the linear regression of abundance (predicted vs. actual)
# Function is useful since it is able to use any of the taxonomic hierarchy in 
# the taxonomy file for further study
regressionScatter<-function(outBacteria, bacteria, name){
  df<-data.frame("Predicted"=predict(outBacteria), "Actual"=bacteria)
  graph<-ggplot(df, 
                aes(x=Predicted, 
                    y=Actual))+
                geom_point()+
                xlab(glue("Predicted {name} abundance"))+
                ylab(glue("Actual {name} abundance"))+
                geom_abline(slope=1, intercept=0)+
                theme_classic(base_size=20)+
                stat_smooth(method="lm")
  print(graph)
  ggsave(glue("Figures/predVsAct{name}Abundance.png"))
}

# Data Import -------------------------------------------------------------

fileNames<-c("about", "from16S", "fromITS", "taxonomy16S", "taxonomyITS")
fileLocation<-c("None", "None", "None", "None", "None") # Removed.

for (i in 1:length(fileNames)){
  loadDataToTxt(fileNames[i], fileLocation[i])
}

#Extract ASV from repective data frames and vectorise it.
vec<-as.vector(from16S$ASV)
taxonomy16S<-taxonomy16S %>% filter(ASV %in% vec)

vec<-as.vector(fromITS$ASV)
taxonomyITS<-taxonomyITS %>% filter(ASV %in% vec)


# py2RNetworkX package requires an input of a data frame as a .txt file.
phylum16S<-groupBy(from16S[,-dim(from16S)[2]], taxonomy16S, "phylum")
write.table(phylum16S, "Data/phylum16S.txt") 

phylumITS<-groupBy(fromITS[,-dim(fromITS)[2]], taxonomyITS, "phylum")
write.table(phylumITS, "Data/phylumITS.txt")

# Creates a copy of phylum16S df and removes specific column names from the df 
# for another four graphs
reducedPhylum16S<-phylum16S
reducedPhylum16S<-reducedPhylum16S[ ,
                      -which(colnames(phylum16S) %in% c("Actinobacteriota",
                        "Bacteroidota", "Proteobacteria", "Myxococcota",
                        "Verrucomicrobiota", "Acidobacteriota"))]
write.table(reducedPhylum16S, "Data/reducedPhylum16S.txt")

reducedPhylumITS<-phylumITS
reducedPhylumITS<-reducedPhylumITS[,
                      -which(colnames(phylumITS) %in% c("Ascomycota",
                        "Rozellomycota", "Glomeromycota"))]
write.table(reducedPhylumITS, "Data/reducedPhylumITS.txt")


# Data Analysis and Graphs ------------------------------------------------

# Knowledge graph (KG) was used due to the spatial power it holds, it allows the 
# user to weight the positions of the bacteria closer to the sites they are 
# found at most frequently, which allows an intuitive view of the spatial 
# features of the community, essentially creating a community ordination
createGraph(df_location="Data/phylum16S.txt", 
        phyla=colnames(phylum16S),
        about_location="Data/Tidy_sampleData_Norway_root_endophytes_York.txt",
        save_location="Figures/knowledgeGraphBacteria.png") 

createGraph(df_location="Data/phylumITS.txt", 
        phyla=colnames(phylumITS),
        about_location="Data/Tidy_sampleData_Norway_root_endophytes_York.txt",
        save_location="Figures/knowledgeGraphFungi.png")


# Reduced Phylum dataframes include only the selected phyla, ones outside the
# middle cluster
createGraph(df_location="Data/reducedPhylum16S.txt", 
          phyla=colnames(reducedPhylum16S),
          about_location="Data/Tidy_sampleData_Norway_root_endophytes_York.txt",
          save_location="Figures/knowledgeGraphBacteriaReduced.png") 

createGraph(df_location="Data/reducedPhylumITS.txt", 
          phyla = colnames(reducedPhylumITS),
          about_location="Data/Tidy_sampleData_Norway_root_endophytes_York.txt",
          save_location="Figures/knowledgeGraphFungiReduced.png")


# Make sure observations of data for ANOVA follow the three assumptions:

# 1. Observations are independent by assumption, data samples were taken far from
# each other there species populations would not interact. The experiment was 
# done with a factoral design, hence observations are independent.

# 2. Observations and normally distributed in every sample verified using a 
# Quartile-Quartile (QQ) plot
qqnorm(phylum16S, pch=1, frame=FALSE)
qqline(phylum16S, col="steelblue", lwd=2)
# The Q-Q plot shows that actually our results aren't normally distributed, 
# however analysis using ANOVA will go ahead anyway.

# 3. The variances is checked by plotting the anova. It fits the assumption
# as the red line is horizontal and doesn't curve much.
plot(res.aov, which=1)

# Call function to plot a boxplot and print two way ANOVA results, TukeyHSD
# results on both the bacteria and fungi phyla 
boxPlotTwoWayAnova(phylumDF=phylum16S)
ggsave("Figures/boxPlotBacteriaShannon.png")
boxPlotTwoWayAnova(phylumDF=phylumITS)
ggsave("Figures/boxPlotFungiShannon.png")


# Group Bacteria and Fungi by class to narrow the results returned
class16S<-groupBy(from16S[,-dim(from16S)[2]], taxonomy16S, "class")
classITS<-groupBy(fromITS[,-dim(fromITS)[2]], taxonomyITS, "class")

# ids and idnums used to call the specific name and row from other data frames
ids<-rownames(class16S)
idNums<-which(about$site %in% ids)

# Calls function to find p values of the bacterial class data frame
pvalsViaRegression(class16S)

# Creates a data frame with the input variables of the linear model 
# (precipitation, temperature) and the actual output values to be plotted against
# the linear regression's predicted values (abundance)
df<-data.frame(bacilli=class16S[, "Bacilli"],
               alpha=class16S[, "Alphaproteobacteria"], 
               precip=about[idNums,"P.grid"],
               temp=about[idNums, "T.grid"])


# Finds p values of Bacilli using Linear Regression, input variables of
# precipitation and temperature, output is abundance
outBacilli<-lm(bacilli ~ P.grid + T.grid, data=df)
outAlpha<-lm(alpha ~ P.grid + T.grid, data=df)
summary(outBacilli)
summary(outAlpha)


# Calls function to plot actual vs. predicted for Bacilli and Alphaproteobacteria
regressionScatter(outBacilli, bacilli, "Bacilli")
regressionScatter(outAlpha, alpha, "Alpha")

# Creating a data frame which includes the order of bacteria as well as the 
# temperature and precipitation
alphaOnly<-from16S[which(taxonomy16S[rownames(from16S), 
                  "class"]=="Alphaproteobacteria"),-dim(from16S)[2]]
order16SAlpha<-groupBy(alphaOnly, taxonomy16S, "order")


threshold<-dim(order16SAlpha)[1]
order16SAlpha<-order16SAlpha[,which(colSums(order16SAlpha)>threshold)]

df<-data.frame(order16SAlpha, 
               precip=about[idNums,"P.grid"],
               temp=about[idNums, "T.grid"])

# Finds p values of Rhizobiales and Sphingomonadales using a Linear Regression as
# before with Bacilli and Alphaproteobacteria

outRhizobiales=lm(Rhizobiales ~ P.grid + T.grid, data=df)
outSphingomonadales=lm(Sphingomonadales ~ P.grid + T.grid, data=df)
summary(outRhizobiales)
summary(outSphingomonadales)

# Calls function to plot actual vs. predicted for Rhizobiales and 
# Sphingomonadales
regressionScatter(outRhizobiales, df[,"Rhizobiales"] ,"Rhizobiales")
regressionScatter(outSphingomonadales, df[,"Sphingomonadales"],
                  "Sphingomonadales")


# Check for statistically significant relationships between fungi frequency and
# temperature precipitation 
ids<-rownames(classITS)
idNums<-which(about$site %in% ids)

# Find p values using a linear regression model for the abundance of the classes 
# of the fungi, with precipitation and temperature as the input variables and 
# abundance as the out

pvalsViaRegression(classITS)
# All p values are above 0.05. Therefore no significant relationships to graph.

