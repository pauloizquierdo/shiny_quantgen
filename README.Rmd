---
title: "Shyny app for quantitative genetics"
#author: "Paulo"
#date: "4/11/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
_________________________________

Created by Sikta Das Adhikari, Claudia Miranda, Joanne Thompson and Paulo Izquierdo, for the a Quantitative Genetics workshop website.

### Dependencies
_________________________________

R and RStudio

Additionally, make sure tidyverse, shiny, rrBLUP and qqman packages are installed.

>[Quantitative genetics workshop videos](https://youtube.com/playlist?list=PLOb4571zCOd8rnWQOTMGnSx5bncpGr9W6)

>[Quantitative genetics workshop website](https://pauloizquierdo.github.io/Quantitative_Genetics/)

_________________________________

```{r echo=FALSE}

hist(rnorm(100000), breaks = 100, 
     col="salmon", xlab = "", main = "")



```


