#' Pick features for cell type prediction with caret by LOOCV
#'
#' @import caret
#' @importFrom rlang .data
#'
#' @param betas A beta matrix of reference DNA methylation data that will be used to select for features for cell type prediction
#' @param meta A data frame of phenotype data, specifying which of the reference samples are of what cell types
#' @param ct A vector of characters specifying the cell types to deconvolute
#' @param caretMods A vector of characters, selecting the models to use to pick the cell type prediction features, options include "lasso", "EL", "BLR", "CART", "RF", "GBM", "PLDA", "GAnRF", "GAnNB", "GAnSVM", and "GAnNN"
#' @param filterK An integer, representing the number probes to input to the machine learning algorithms with top T test probes
#' @param seed An integer specifying the seed for reproducibility
#' @param plot A Boolean specifying whether to plot a heatmap showing the clustering performance using the probes selected
#' @param verbose A Boolean specifying whether the function should be verbose or not

pickCompProbesCaretLOOCV <- function(betas, meta, ct,
                                     caretMods = c("lasso", "EL", "BLR", "CART", "RF", "GBM", "PLDA", "GAnRF", "GAnNB", "GAnSVM", "GAnNN"),
                                     filterK = 10000, seed = 1234, plot = TRUE, verbose = TRUE) {
    df <- as.matrix(betas)
    pd <- as.data.frame(meta)

    ## only keep the cell type you're estimating
    keep <- which(pd$cellType %in% ct)
    pd <- pd[keep, ]
    df <- df[, keep]
    pd$cellType <- factor(pd$cellType, levels = ct) # make cell type a factor

    ## Find the probes with significant different methylation status across cell types
    tIndexes <- splitit(pd$cellType)
    tIndexes <- tIndexes[lapply(tIndexes, length) > 0]

    ## set seeds for reproducibility
    set.seed(seed)
    seeds5 <- vector(mode = "list", length = ncol(df)+1) # LOOCV - change if the sampling method is altered in trainControl below
    for(i in 1:ncol(df)) seeds5[[i]] <- sample.int(1000, 5)
    seeds5[[ncol(df)+1]] <- sample.int(1000, 1)
    seeds9 <- vector(mode = "list", length = ncol(df)+1) # LOOCV - change if the sampling method is altered in trainControl below
    for(i in 1:ncol(df)) seeds9[[i]] <- sample.int(1000, 9)
    seeds9[[ncol(df)+1]] <- sample.int(1000, 1)
    seeds10 <- vector(mode = "list", length = ncol(df)+1) # LOOCV - change if the sampling method is altered in trainControl below
    for(i in 1:ncol(df)) seeds10[[i]] <- sample.int(1000, 10)
    seeds10[[ncol(df)+1]] <- sample.int(1000, 1)
    seeds400 <- vector(mode = "list", length = ncol(df)+1) # LOOCV - change if the sampling method is altered in trainControl below
    for(i in 1:ncol(df)) seeds400[[i]] <- sample.int(1000, 400)
    seeds400[[ncol(df)+1]] <- sample.int(1000, 1)

    seedsGA <- sample.int(1000, ncol(df)+1)
    seedsGA2 <- vector(mode = "list", length = ncol(df)+1) # LOOCV - change if the sampling method is altered in trainControl below
    for(i in 1:ncol(df)) seedsGA2[[i]] <- sample.int(1000, 3)
    seedsGA2[[ncol(df)+1]] <- sample.int(1000, 1)

    control5 <- caret::trainControl(method = "LOOCV",
                                    number = ncol(df),
                                    verboseIter = F,
                                    seeds = seeds5)
    control10 <- caret::trainControl(method = "LOOCV",
                                     number = ncol(df),
                                     verboseIter = F,
                                     seeds = seeds10)
    control9 <- caret::trainControl(method = "LOOCV",
                                    number = ncol(df),
                                    verboseIter = F,
                                    seeds = seeds9)
    control400 <- caret::trainControl(method = "LOOCV",
                                      number = ncol(df),
                                      verboseIter = F,
                                      seeds = seeds400)
    controlGA <- caret::trainControl(method = "LOOCV",
                                     number = ncol(df),
                                     verboseIter = F,
                                     seeds = seedsGA2)

    ## Setting plotting colors
    if (plot){
        pltct <- c(Bcell = "#980043FF", Bmem = "#CE1256FF", Bnv = "#DF65B0FF",
                   CD4T = "#662506FF", CD4mem = "#993404FF", CD4nv = "#FE9929FF",
                   CD8T = "#B30000FF", CD8mem = "#CB181DFF", CD8nv = "#FB6A4AFF",
                   Treg = "#FF847CFF", Mono = "#2F8AC4", NK = "#764E9F",
                   Gran = "#99C945", Neu = "#52BCA3", Eos = "#24796C",
                   Bas = "#A8DDB5FF", nRBC = "#67000DFF", PBMC = "#A5AA99", WBC = "#252525FF")
        pltct <- pltct[names(pltct) %in% ct]
        anncolors <- list(cellType = pltct)
    }

    ## Pick probes based on specified machine learning methods
    probeList <- lapply(ct, function(ctType) {
        x <- tIndexes[[ctType]]
        ctIndex <- rep("N", ncol(df))
        ctIndex[x] <- "Y"
        ctIndex <- as.factor(ctIndex)
        tout <- genefilter::rowttests(df, ctIndex)
        tout.top <- rownames(tout)[order(tout$p.value, decreasing = F)][1:filterK] # pick filterK number of probes as the first pass
        df <- df[tout.top, ] %>% t()

        out <- vector(mode = 'list', length = length(caretMods) + 1)
        names(out) <- c(caretMods, "tTestTopK")

        # lasso:
        # Use L1 regularization to optimize shrinkage of coefficients uncorrelated with outcome of interest
        if ("lasso" %in% caretMods){
            requireNamespace("glmnet")
            if(verbose) cat(paste0("Running lasso for feature selection of ", ctType, ".\n"))
            lassotune <- expand.grid(alpha = 1,
                                     lambda = seq(0.1, 1, length = 10))
            lassoout <- caret::train(x = df,
                                     y = ctIndex,
                                     method = "glmnet",
                                     tuneGrid = lassotune,
                                     trControl = control10)
            if(verbose) cat(paste0("Number of features selected by lasso for ", ctType, ": ", length(caret::predictors(lassoout)), "\n"))
            if (plot & length(caret::predictors(lassoout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(lassoout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("lasso - ", ctType))
            }
            out[["lasso"]] <- list(coefs = caret::predictors(lassoout), performance = lassoout$resample)
        }
        # Elastic net:
        # Elastic Net is a linear regression-based method that can perform feature selection by using L1 and L2 regularization.
        # The L1 regularization can shrink coefficients to zero and eliminate irrelevant features, while the L2 regularization can help to reduce overfitting.
        if ("EL" %in% caretMods){
            requireNamespace("glmnet")
            if(verbose) cat(paste0("Running elastic net (EL) for feature selection of ", ctType, ".\n"))
            ELtune <- expand.grid(alpha = seq(0, 1, length = 10),
                                  lambda = seq(0.1, 1, length = 10))
            ELout <- caret::train(x = df,
                                  y = ctIndex,
                                  method = "glmnet",
                                  tuneGrid = ELtune,
                                  trControl = control10)
            if(verbose) cat(paste0("Number of features selected by EL for ", ctType, ": ", length(caret::predictors(ELout)), "\n"))
            if (plot & length(caret::predictors(ELout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(ELout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("elastic net - ", ctType))
            }
            out[["EL"]] <- list(coefs = caret::predictors(ELout), performance = ELout$resample)
        }
        # Random Forest:
        # Random Forest is a decision tree-based ensemble method that can perform feature selection by evaluating the importance of each feature in the model.
        # The Gini importance or Mean Decrease Impurity (MDI) are two ways to calculate feature importance in Random Forest.
        if ("RF" %in% caretMods){
            if(verbose) cat(paste0("Running random forest (RF) for feature selection of ", ctType, ".\n"))
            RFtune <- expand.grid(mtry = seq(10, 1500, length = 10))
            RFout <- caret::train(x = df,
                                  y = ctIndex,
                                  method = "rf",
                                  ntree = 500,
                                  nodesize = 5,
                                  tuneGrid = RFtune,
                                  trControl = control10)
            if(verbose) cat(paste0("Number of features selected by RF for ", ctType, ": ", length(caret::predictors(RFout)), "\n"))
            if (plot & length(caret::predictors(RFout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(RFout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("random forest - ", ctType))
            }
            out[["RF"]] <- list(coefs = caret::predictors(RFout), performance = RFout$resample)
        }
        # Boosted logistic regression:
        # Boosted logistic regression is a machine learning algorithm that combines logistic regression with boosting, a technique for iteratively improving the performance of weak learners.
        # In boosted logistic regression, a weak logistic regression model is trained on the data, and then additional weak logistic regression models are trained to correct the errors of the previous models.
        # Each subsequent model focuses more on the data points that were misclassified by the previous models, thus improving the overall accuracy of the model.
        if ("BLR" %in% caretMods){
            requireNamespace("caTools")
            if(verbose) cat(paste0("Running boosted logistic regression (BLR) for feature selection of ", ctType, ".\n"))
            BLRout <- caret::train(x = df,
                                   y = ctIndex,
                                   method = "LogitBoost",
                                   tuneLength = 5,
                                   trControl = control5)
            if(verbose) cat(paste0("Number of features selected by BLR for ", ctType, ": ", length(caret::predictors(BLRout)), "\n"))
            if (plot & length(caret::predictors(BLRout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(BLRout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("boosted logistic regression - ", ctType))
            }
            out[["BLR"]] <- list(coefs = caret::predictors(BLRout), performance = BLRout$resample)
        }
        # Classification And Regression Tree:
        # Decision Trees can perform feature selection by evaluating the information gain or Gini impurity of each feature when deciding how to split the data.
        # The important features are those that lead to the highest reduction in impurity or the highest information gain.
        if ("CART" %in% caretMods){
            requireNamespace("rpart")
            if(verbose) cat(paste0("Running Classification And Regression Tree (CART) for feature selection of ", ctType, ".\n"))
            CARTout <- caret::train(x = df,
                                    y = ctIndex,
                                    method = "rpart",
                                    tuneLength = 5,
                                    trControl = control5)
            if(verbose) cat(paste0("Number of features selected by CART for ", ctType, ": ", length(caret::predictors(CARTout)), "\n"))
            if (plot & length(caret::predictors(CARTout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(CARTout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("classification and regression tree - ", ctType))
            }
            out[["CART"]] <- list(coefs = caret::predictors(CARTout), performance = CARTout$resample)
        }
        # Gradient Boosting Machines:
        # GBM is another ensemble-based method that can perform feature selection by iteratively building trees that focus on predicting the errors of previous trees.
        # In this process, GBM can identify important features that improve the model's performance.
        if ("GBM" %in% caretMods){
            requireNamespace("xgboost")
            if(verbose) cat(paste0("Running Gradient Boosting Machines (GBM) for feature selection of ", ctType, ".\n"))
            GBMout <- caret::train(x = as.data.frame(df),
                                   y = ctIndex,
                                   method = "xgbDART", # xgbDART
                                   metric = "Accuracy",
                                   trControl = control400)
            if(verbose) cat(paste0("Number of features selected by GBM for ", ctType, ": ", length(caret::predictors(GBMout)), "\n"))
            if (plot & length(caret::predictors(GBMout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(GBMout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("gradient boosting machines - ", ctType))
            }
            out[["GBM"]] <- list(coefs = caret::predictors(GBMout), performance = GBMout$resample)
        }
        # Penalized Linear Discriminant Analysis (PLDA):
        # LDA is a technique for dimensionality reduction, meaning that it reduces the number of input features to a smaller set of new features that are most discriminative for the classification task.
        # Specifically, LDA aims to find a linear combination of the input features that maximizes the separation between the classes, while minimizing the variance within each class.
        if ("PLDA" %in% caretMods){
            requireNamespace("rqPen")
            if(verbose) cat(paste0("Running Penalized Linear Discriminant Analysis (PLDA) for feature selection of ", ctType, ".\n"))
            set.seed(seed)
            PLDAout <- train(x = df,
                             y = ctIndex,
                             method = "PenalizedLDA",
                             tuneLength = 5,
                             trControl = control5)
            if(verbose) cat(paste0("Number of features selected by PLDA for ", ctType, ": ", length(caret::predictors(PLDAout)), "\n"))
            if (plot & length(caret::predictors(PLDAout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(PLDAout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("genetic algorithm + linear discriminant analysis - ", ctType))
            }
            out[["PLDA"]] <- list(coefs = caret::predictors(PLDAout), performance = PLDAout$resample)
        }
        # Simulated annealing + Neural Networks:
        # Neural Networks can perform feature selection by learning a set of weights that indicate the importance of each feature in predicting the output.
        # The weights can be used to identify the most important features, and some neural network architectures, such as Sparse Autoencoders, are designed to learn sparse representations, which effectively perform feature selection.
        if ("GAnNN" %in% caretMods){
            requireNamespace("nnet")
            if(verbose) cat(paste0("Running Genentic Annealing (GA) + Neural Networks (NN) for feature selection of ", ctType, ".\n"))
            # NNout <- train(x = df,
            #                y = ctIndex,
            #                method = "nnet",
            #                tuneLength = 3,
            #                trControl = control9)
            GAnNNout <- caret::gafs(x = df,
                                    y = ctIndex,
                                    iters = 1,
                                    gafsControl = gafsControl(method = "LOOCV",
                                                              functions = caretGA,
                                                              seeds = seedsGA,
                                                              verbose = TRUE),
                                    method = "nnet",
                                    tuneLength = 3,
                                    trControl = control9)
            if(verbose) cat(paste0("Number of features selected by GAnNN for ", ctType, ": ", length(caret::predictors(GAnNNout)), "\n"))
            if (plot & length(caret::predictors(GAnNNout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(GAnNNout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("genetic algorithm + neural networks - ", ctType))
            }
            out[["GAnNN"]] <- list(coefs = caret::predictors(GAnNNout), performance = GAnNNout$resample)
        }
        # Simulated annealing + Naive Bayes:
        # Naive Bayes is a probabilistic classification method that can perform feature selection by calculating the probability of each feature given the class.
        # The most important features are those that have the highest probability of occurring in the positive class.
        if ("GAnNB" %in% caretMods){
            requireNamespace("naivebayes")
            if(verbose) cat(paste0("Running Genentic Annealing (GA) + Naive Bayes (NB) for feature selection of ", ctType, ".\n"))
            # NBout <- train(x = df,
            #                y = ctIndex,
            #                method = "nb",
            #                tuneLength = 5,
            #                trControl = control5)
            GAnNBout <- caret::gafs(x = df,
                                    y = ctIndex,
                                    iters = 1,
                                    gafsControl = gafsControl(method = "LOOCV",
                                                              functions = caretGA,
                                                              seeds = seedsGA),
                                    method = "nb",
                                    tuneLength = 3,
                                    trControl = controlGA)
            if(verbose) cat(paste0("Number of features selected by GAnNB for ", ctType, ": ", length(caret::predictors(GAnNBout)), "\n"))
            if (plot & length(caret::predictors(GAnNBout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(GAnNBout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("ggenetic algorithm + naive Bayes - ", ctType))
            }
            out[["GAnNB"]] <- list(coefs = caret::predictors(GAnNBout), performance = GAnNBout$resample)
        }
        # Support Vector Machines (SVM):
        # SVM is a classification method that can perform feature selection by maximizing the margin between two classes.
        # In this process, SVM can identify important features that contribute the most to the classification task.
        if ("GAnSVM" %in% caretMods){
            requireNamespace("kernlab")
            if(verbose) cat(paste0("Running Genentic Annealing (GA) + Support Vector Machines (SVM) for feature selection of ", ctType, ".\n"))
            # SVMout <- train(x = df,
            #                 y = ctIndex,
            #                 method = "svmLinear",
            #                 tuneLength = 5,
            #                 trControl = control5)
            GAnSVMout <- caret::gafs(x = df,
                                     y = ctIndex,
                                     iters = 1,
                                     gafsControl = gafsControl(method = "LOOCV",
                                                               functions = caretGA,
                                                               seeds = seedsGA),
                                     method = "svmLinear",
                                     tuneLength = 3,
                                     trControl = controlGA)
            if(verbose) cat(paste0("Number of features selected by GAnSVM for ", ctType, ": ", length(caret::predictors(GAnSVMout)), "\n"))
            if (plot & length(caret::predictors(GAnSVMout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(GAnSVMout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("genetic algorithm + support vector machines - ", ctType))
            }
            out[["GAnSVM"]] <- list(coefs = caret::predictors(GAnSVMout), performance = GAnSVMout$resample)
        }
        # Simulated annealing + Random forest
        if ("GAnRF" %in% caretMods){
            if(verbose) cat(paste0("Running Genentic Annealing (GA) + Random Forest (RF) for feature selection of ", ctType, ".\n"))
            set.seed(seed)
            GAnRFout <- caret::gafs(x = df,
                                    y = ctIndex,
                                    iters = 1,
                                    gafsControl = gafsControl(method = "LOOCV",
                                                              functions = caretGA,
                                                              seeds = seedsGA),
                                    method = "rf",
                                    tuneLength = 3,
                                    trControl = controlGA)
            if(verbose) cat(paste0("Number of features selected by GAnRF for ", ctType, ": ", length(caret::predictors(GAnRFout)), "\n"))
            if (plot & length(caret::predictors(GAnRFout)) > 1){
                plt <- t(df)
                pheatmap::pheatmap(plt[rownames(plt) %in% caret::predictors(GAnRFout), ],
                                   show_colnames = F,
                                   show_rownames = F,
                                   annotation_col = pd[, "cellType", drop = F],
                                   annotation_colors = anncolors,
                                   main = paste0("genetic algorithm + random forest - ", ctType))
            }
            out[["GAnRF"]] <- list(coefs = caret::predictors(GAnRFout), performance = GAnRFout$resample)
        }
        # Compare performance
        # tests <- ls()[ls() %in% paste0(caretMods, "out")]
        # allTests <- lapply(1:length(tests), function(x){
        #     return(get(tests[x]))
        # })
        # names(allTests) <- tests
        # resamps <- resamples(allTests)
        # out[["summary"]] <- resamps
        # print(summary(resamps))
        # GAout <- gafs(x = df,
        #             y = as.factor(ctIndex),
        #             iters = 50)
        # smallest <- pickSizeBest(example, metric = "RMSE", maximize = FALSE)
        out[["tTestTopK"]] <- tout.top
        return(out) # For each cell type, compare the mean DNAm level in this cell type vs all others
    })
    names(probeList) <- ct

    combinedprobes <- lapply(caretMods, function(fs){ # Combine all probes
        probes <- lapply(probeList, function(x){
            out <- x[[fs]][["coefs"]]
            return(out)
        })
        trainingProbes <- unlist(probes)
        return(trainingProbes)
    }) %>% do.call(c, .) %>% unique()
    dfComb <- df[combinedprobes, ] # Subset down the reference data to the selected probes - remove duplicated CpG since they generate the same coefficients
    mod <- stats::model.matrix(~ pd$cellType - 1) %>% as.data.frame()
    colnames(mod) <- levels(pd$cellType)
    form <- stats::as.formula(sprintf("bt ~ %s - 1", paste(ct, collapse = " + ")))
    coefComb <- sapply(1:nrow(dfComb), function(probe){
        mod$bt <- df[probe, ]
        fit <- stats::lm(form, data = mod[stats::complete.cases(mod), ]) # Remove samples with missing methylation values
        fitCoef <- fit$coef
        return(fitCoef)
    }) %>% t() %>% as.data.frame()
    rownames(coefComb) <- rownames(dfComb)

    probeCoefs <- lapply(caretMods, function(fs){
        probes <- lapply(probeList, function(x){
            out <- x[[fs]][["coefs"]]
            return(out)
        })
        trainingProbes <- unlist(probes)
        ## Call a linear model with the selected probes and calculate the weights of each cell type
        df <- df[unique(trainingProbes), ] # Subset down the reference data to the selected probes - remove duplicated CpG since they generate the same coefficients
        mod <- stats::model.matrix(~ pd$cellType - 1) %>% as.data.frame()
        colnames(mod) <- levels(pd$cellType)
        form <- stats::as.formula(sprintf("bt ~ %s - 1", paste(ct, collapse = " + ")))
        coef <- sapply(1:nrow(df), function(probe){
            mod$bt <- df[probe, ]
            fit <- stats::lm(form, data = mod[stats::complete.cases(mod), ]) # Remove samples with missing methylation values
            fitCoef <- fit$coef
            return(fitCoef)
        }) %>% t() %>% as.data.frame()
        rownames(coef) <- rownames(df)
        return(coef)
    })
    names(probeCoefs) <- caretMods
    probeCoefs[["combined"]] <- coefComb

    return(list(probeCoefs = probeCoefs, probeList = probeList))
}
