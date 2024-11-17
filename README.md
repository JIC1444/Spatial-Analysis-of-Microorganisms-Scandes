# Project Overview
This analysis was done as part of a larger project, data was provided from researchers from the Scandes Mountains in Norway on the bacteria and fungi present at 12 different sites in the mountains, analyse the given data and contribute a brand new insight, filling a gap in the literature.

This project persued the question "Is there a correlation between the receding tree line* and the microbes found there?". This follows from the master's thesis [1] and shows through the choice of the figures in the paper - with spatial elements, trying to determine if there are a group (or lack thereof) of bacteria and or fungi which contribute to the increasing bare space on the Scandes mountains.

_* A tree line is the ring of trees adjacent to the bare top of a mountain - it is a metric used to track the progress of global warming._

###### _Data Privacy: Due to limitations with the source of the data, it cannot be released in conjunction with the analysis file, however my figures and code can be shared._

# Figures

## Figure 1: Knowledge graph depicting the environments favoured by each phyla of bacteria/fungi
<img src="https://github.com/user-attachments/assets/449b5c32-93a3-4626-857e-cb44150d1574" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/c910eb4b-39b3-4d24-a4f5-50b80366909a" width="50%" height="50%">
<img src="https://github.com/user-attachments/assets/14a489ee-baea-435b-93ae-cc3aae3a629c" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/6bcb54d8-9700-4148-bbae-0e59fa4d3bc1" width="50%" height="50%">

The graph shows three vertical lines, which align with the elevational trend for these sites; the from left to right the sites increase in elevation. The closer a phylum node to a location node, the greater the count found there, the highlighted edges show a phylum's strongest connection to a site. The top row of graphs is the bacteria phyla, with the left image displaying all of the phyla in the study. The right image has had all of the central phyla removed from the graph, to highlight the 'outliers', and persue scientific literature on these bacteria and their potential effect on the treeline.

The motivation for this figure was to capture all of the data's elements into one. A PCA would have only covered the temperature and precipitation within the data, missing the elevation. The PCA also would have been unstructured whereas this graph has edges between the nodes, making the connection clear between a site and its constituent microbes.

I made a custom R package to essentially port the functionality from NetworkX in Python to make a graph in R. The reasoning behind this was that if I wanted to do any serious machine or deep learning on this dataset, then having this NetworkX object allows very easy conversion back to Python (where the model would be defined) then is able to be fed in to said model right away. I however, decided that linear or random forest regression was more suited to the data at hand (Figure 3).

```{r}
# Knowledge graph (KG) was used due to the spatial power it holds, it allows the 
# user to weight the positions of the bacteria closer to the sites they are 
# found at most frequently, which allows an intuitive view of the spatial 
# features of the community, essentially creating a community ordination.
createGraph(df_location="Data/phylum16S.txt", 
        phyla=colnames(phylum16S),
        about_location="Data/Tidy_sampleData_re.txt",
        save_location="Figures/knowledgeGraphBacteria.png") 

createGraph(df_location="Data/phylumITS.txt", 
        phyla=colnames(phylumITS),
        about_location="Data/Tidy_sampleData_re.txt",
        save_location="Figures/knowledgeGraphFungi.png")


# Reduced Phylum dataframes include only the selected phyla, ones outside the
# middle cluster
createGraph(df_location="Data/reducedPhylum16S.txt", 
          phyla=colnames(reducedPhylum16S),
          about_location="Data/Tidy_sampleData_re.txt",
          save_location="Figures/knowledgeGraphBacteriaReduced.png") 

createGraph(df_location="Data/reducedPhylumITS.txt", 
          phyla = colnames(reducedPhylumITS),
          about_location="Data/Tidy_sampleData_re.txt",
          save_location="Figures/knowledgeGraphFungiReduced.png")

```

The Python code for the 'createGraph' function can be found in the appendix.

## Figure 2: Boxplot of shannon index vs. the site conditions
<img src="https://github.com/user-attachments/assets/be893bef-8c29-493a-b887-b0410ebbcf5c" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/f105f6e6-8caf-456d-85d0-8c34381631c9" width="50%" height="50%">

These graphs show the previous figure's information in a more statistical manner, with mean, median and range. This highlights whether there is any significant statistical relationship between the Shannon diversity index to the temperature and the precipitiation of the sites. After ANOVA testing, there were a few significant relationships, these, I hoped would reveal more 
these could then be evaluated with further testing within the data.

First of all, an ANOVA test has three requirements for the data.
```{r} 
# Make sure observations of data for ANOVA follow the three assumptions:

# 1. Observations are independent by assumption, data samples were taken far from
# each other there species populations would not interact. The experiment was 
# done with a factoral design, hence observations are independent.

# 2. Observations and normally distributed in every sample verified using a 
# Quartile-Quartile (QQ) plot
qqnorm(phylum16S, pch=1, frame=FALSE)
qqline(phylum16S, col="steelblue", lwd=2)
# The Q-Q plot shows that actually our data aren't normally distributed

# 3. The variance is checked by plotting the anova. It fits the assumption
# as the red line is horizontal and doesn't curve much.
plot(res.aov, which=1)
```
The data fails the second test, however we ignore this and proceed as the Q-Q plot only showed a minor deviation from the normal distribution.

```{r}
# Plot a boxplot and print two way ANOVA results, TukeyHSD
# results on both the bacteria and fungi phyla 
boxPlotTwoWayAnova(phylumDF=phylum16S)
ggsave("Figures/boxPlotBacteriaShannon.png")
boxPlotTwoWayAnova(phylumDF=phylumITS)
ggsave("Figures/boxPlotFungiShannon.png")
```

## Figure 3: Acutal vs. predicted bacteria abundance
<img src="https://github.com/user-attachments/assets/db4a3997-8f3a-47ab-9f76-ff4cc93ad671" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/75f288fa-2bc2-4bb3-9012-b65775d9f7d0" width="50%" height="50%">
<img src="https://github.com/user-attachments/assets/79091112-0cd5-4efb-b038-e64c8a926ef1" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/0cf4dee2-72c7-4624-bd2a-a3c6e242a114" width="50%" height="50%">

This graph shows whether there actually is a relationship between the temperature and precipitation of a site and the microbes found at it - these four bacteria were selected due to their high statistical significance. These showed promise toward the conclusion in [1], where it was stated that there is a relationship between the temperature and the precipitation of the site, to whether it was above, below or on the tree line. This then may coincide to certain types of micro-organisms which are or are not found in that specific temperature-precipitation combination.

This code probes the data to find any significant relationships in the higher orders of the bacteria.
```{r}
# Calls function to find p values of the bacterial class data frame
pvalsViaRegression(class16S)
# Bacilli and Alphaproteobacteria show significance.

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

```

This code then looks for any significant relationships in the fungal orders.
```{r}
# Check for statistically significant relationships between fungi frequency and
# temperature precipitation 
ids<-rownames(classITS)
idNums<-which(about$site %in% ids)

# Find p values using a linear regression model for the abundance of the classes 
# of the fungi, with precipitation and temperature as the input variables and 
# abundance as the out

pvalsViaRegression(classITS)
# All p values are above 0.05. Therefore no significant relationships to graph.
```

No significant relationships found in the fungi! Hence why there are bacteria only in Figure 3.

### Conclusion
Since this project spanned 2 months and used research from others it cannot be said for certain whether these microbes contribute to the tree line in Norway, however there were some data which pointed toward specific bacterial communities contributing to the growth of young saplings. There are grounds for further study into this theory of the relationship between the tree line and the temperature, precipitation and the microbes found there.


## References
[1] Tingstad, L. (2012) Tree seedling establishment along climatic gradients from lowland to alpine areas in western Norway, Department of Nature and Resource Management, Master’s Thesis


## Appendix: R functions
```{r}
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
```
