#load (and install, if necessary) packages
tryCatch(library("ggplot2"), 
         error = function(e){
           install.packages(pkgs =  "ggplot2", 
                            repos = 'http://cran.us.r-project.org')
           library("ggplot2")
         })
tryCatch(library("caret"), 
         error = function(e){
           install.packages(pkgs =  "caret", 
                            repos = 'http://cran.us.r-project.org')
           library("caret")
         })
tryCatch(library("openxlsx"), 
         error = function(e){
           install.packages(pkgs =  "openxlsx", 
                            repos = 'http://cran.us.r-project.org')
           library("openxlsx")
         })
tryCatch(library("doParallel"), 
         error = function(e){
           install.packages(pkgs =  "doParallel", 
                            repos = 'http://cran.us.r-project.org')
           library("doParallel")
         })
tryCatch(library("rms"), 
         error = function(e){
           install.packages(pkgs =  "rms", 
                            repos = 'http://cran.us.r-project.org')
           library("rms")
         })
tryCatch(library("dplyr"), 
         error = function(e){
           install.packages(pkgs =  "dplyr", 
                            repos = 'http://cran.us.r-project.org')
           library("dplyr")
         })
tryCatch(library("survival"), 
         error = function(e){
           install.packages(pkgs =  "survival", 
                            repos = 'http://cran.us.r-project.org')
           library("survival")
         })
tryCatch(library("survminer"), 
         error = function(e){
           install.packages(pkgs =  "survminer", 
                            repos = 'http://cran.us.r-project.org')
           library("survminer")
         })

tryCatch(library("ggthemes"), 
         error = function(e){
           install.packages(pkgs =  "ggthemes", 
                            repos = 'http://cran.us.r-project.org')
           library("ggthemes")
         })

tryCatch(library("survAUC"), 
         error = function(e){
           install.packages(pkgs =  "survAUC", 
                            repos = 'http://cran.us.r-project.org')
           library("survAUC")
         })

tryCatch(library("glmnet"), 
         error = function(e){
           install.packages(pkgs =  "glmnet", 
                            repos = 'http://cran.us.r-project.org')
           library("glmnet")
         })
tryCatch(library("reshape2"), 
         error = function(e){
           install.packages(pkgs =  "reshape2", 
                            repos = 'http://cran.us.r-project.org')
           library("reshape2")
         })

# color, vizualisation
tryCatch(library("scales"), 
         error = function(e){
           install.packages(pkgs =  "scales", 
                            repos = 'http://cran.us.r-project.org')
           library("scales")
         })

tryCatch(library("biomaRt"), 
         error = function(e){
           if (!requireNamespace("BiocManager",
                                 quietly = TRUE))
             install.packages("BiocManager", 
                              repos = 'http://cran.us.r-project.org')
           BiocManager::install("biomaRt",
                                update = FALSE,
                                ask = FALSE)
           library("biomaRt")
         })

tryCatch(library("subSeq"), 
         error = function(e){
           if (!requireNamespace("BiocManager",
                                 quietly = TRUE))
             install.packages("BiocManager", 
                              repos = 'http://cran.us.r-project.org')
           BiocManager::install("subSeq",
                                update = FALSE,
                                ask = FALSE)
           library("subSeq")
         })

tryCatch(library("edgeR"), 
         error = function(e){
           if (!requireNamespace("BiocManager",
                                 quietly = TRUE))
             install.packages("BiocManager", 
                              repos = 'http://cran.us.r-project.org')
           BiocManager::install("edgeR",
                                update = FALSE,
                                ask = FALSE)
           library("edgeR")
         })

tryCatch(library("limma"), 
         error = function(e){
           if (!requireNamespace("BiocManager",
                                 quietly = TRUE))
             install.packages("BiocManager", 
                              repos = 'http://cran.us.r-project.org')
           BiocManager::install("limma",
                                update = FALSE,
                                ask = FALSE)
           library("limma")
         })

tryCatch(library("survcomp"), 
         error = function(e){
           if (!requireNamespace("BiocManager",
                                 quietly = TRUE))
             install.packages("BiocManager", 
                              repos = 'http://cran.us.r-project.org')
           BiocManager::install("survcomp",
                                update = FALSE,
                                ask = FALSE)
           library("survcomp")
         })


# for the colors
tryCatch(library("wesanderson"), 
         error = function(e){
           install.packages("wesanderson", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("wesanderson")
           pal <- wes_palette("Zissou1", 100, type = "continuous")
           
         })


tryCatch(library("ranger"), 
         error = function(e){
           install.packages("ranger", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("ranger")
         })
tryCatch(library("tuneRanger"), 
         error = function(e){
           install.packages("tuneRanger", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("tuneRanger")
         })


tryCatch(library("randomForestSRC"), 
         error = function(e){
           install.packages("randomForestSRC", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("randomForestSRC")
         })

tryCatch(library("ggpubr"), 
         error = function(e){
           install.packages("ggpubr", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("ggpubr")
         })

tryCatch(library("ggrepel"), 
         error = function(e){
           install.packages("ggrepel", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("ggrepel")
         })

tryCatch(library("pdftools"), 
         error = function(e){
           install.packages("pdftools", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("pdftools")
         })

tryCatch(library("VennDiagram"), 
         error = function(e){
           install.packages("VennDiagram", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("VennDiagram")
         })

tryCatch(library("DESeq2"),
         error = function(e){
           if (!requireNamespace("BiocManager",
                                 quietly = TRUE))
             install.packages("BiocManager",
                              repos = 'http://cran.us.r-project.org')
           BiocManager::install("DESeq2",
                                update = FALSE,
                                ask = FALSE)
           library("DESeq2")
         })

tryCatch(library("SIS"), 
         error = function(e){
           install.packages("SIS", repos="http://cran.us.r-project.org", 
                            dependencies=TRUE)
           library("SIS")
         })
