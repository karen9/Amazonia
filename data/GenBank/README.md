__GenBank Data__

All genes were downloaded, aligned and saved in Nexus format using the Pipeline 1 of the R-package [alignTools](https://github.com/oleon12/alignTools). Here is an example:

```
library(dplyr)
library(alignTools)

data <- read.csv("Cebidae.csv", header=T)

#Nexus
multiGenBank(data,TRUE,FALSE)%>%
multiMuscle(write.dna = FALSE)%>%
concatGenes(missing = "?", write.dna = TRUE, write.format = "nexus", filename = "concat")
```
