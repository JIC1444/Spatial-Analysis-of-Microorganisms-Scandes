# Project Overview
This analysis was done as part of a larger project, data was provided from researchers from the Scandes Mountains in Norway on the bacteria and fungi present at 8 different sites in the mountains, analyse the given data and contribute a new insight, filling a gap in the literature.

This project persued the question "Is there a correlation between the receding tree line and the microbes found there?". This follows from the master's thesis [] and shows through the choice of the figures in the paper - with spatial elements, trying to determine if there are a group (or lack thereof) of bacteria and or fungi which contribute to the increasing bare space on the Scandes mountains.

##### _Data Privacy: Due to limitations with the source of the data, it cannot be released in conjunction with the analysis file, however my figures can still be shared on this page._

# Figures

## Figure 1: Knowledge graph depicting the environments favoured by each phyla of bacteria/fungi
<img src="https://github.com/user-attachments/assets/449b5c32-93a3-4626-857e-cb44150d1574" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/c910eb4b-39b3-4d24-a4f5-50b80366909a" width="50%" height="50%">
<img src="https://github.com/user-attachments/assets/14a489ee-baea-435b-93ae-cc3aae3a629c" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/6bcb54d8-9700-4148-bbae-0e59fa4d3bc1" width="50%" height="50%">

The graph shows three vertical lines, which align with the elevational trend for these sites; the from left to right the sites increase in elevation. The closer a phylum node to a location node, the greater the count found there, the highlighted edges show a phylum's strongest connection to a site. The top row of graphs is the bacteria phyla, with the left image displaying all of the phyla in the study. The right image has had all of the central phyla removed from the graph, to highlight the 'outliers', and persue scientific literature on these bacteria and their potential effect on the treeline.

The motivation for this figure was to capture all of the data's elements into one. A PCA would have only covered the temperature and precipitation within the data, missing the elevation. The PCA also would have been unstructured whereas this graph has edges between the nodes, making the connection clear between a site and its constituent microbes.

## Figure 2: Boxplot of shannon index vs. the site conditions
<img src="https://github.com/user-attachments/assets/be893bef-8c29-493a-b887-b0410ebbcf5c" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/f105f6e6-8caf-456d-85d0-8c34381631c9" width="50%" height="50%">

These graphs show the previous figure's information more statistically, with mean, median and range. This highlights whether there is any significant statistical relationship between the Shannon diversity index to the temperature and the precipitiation of the sites. After ANOVA testing, there were a few significant relationships, these could then be evaluated with further testing within the data.


## Figure 3: Acutal vs. predicted bacteria abundance
<img src="https://github.com/user-attachments/assets/db4a3997-8f3a-47ab-9f76-ff4cc93ad671" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/75f288fa-2bc2-4bb3-9012-b65775d9f7d0" width="50%" height="50%">
<img src="https://github.com/user-attachments/assets/79091112-0cd5-4efb-b038-e64c8a926ef1" width="50%" height="50%"><img src="https://github.com/user-attachments/assets/0cf4dee2-72c7-4624-bd2a-a3c6e242a114" width="50%" height="50%">

This graph shows whether there actually is a relationship between the temperature and precipitation of a site and the microbes found at it - these four bacteria were selected due to their high statistical significance. These showed promise toward the conclusion in [1], where it was stated that there is a relationship between the 

### Conclusion
Since this project spanned 2 months and only used research from others it cannot be said for certain whether these microbes contribute to the tree line in Norway, however there are grounds for further study into this theory.


## References
[1]
