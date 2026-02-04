################################################################################

.splitit <- function(x) {
    split(seq(along = x), x)
}

################################################################################

.rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

################################################################################

# By ricardo on Stack Overflow
.loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

################################################################################

### Call the right manifest based on the array type - accessory function
# Example usage:
# manifestv2 <- .getManifest("IlluminaHumanMethylationEPICv2")
.getManifest <- function(outType = c("IlluminaHumanMethylation450k",
                                    "IlluminaHumanMethylationEPIC",
                                    "IlluminaHumanMethylationEPICv2")){
    outType <- match.arg(outType)
    out <- if(requireNamespace(paste0(outType, "manifest"))) {
        library(paste0(outType, "manifest"), character.only = T)
        message(sprintf("manifest package for %s loaded", outType))
        get(paste0(outType, "manifest"))
    } else stop(sprintf("cannot load manifest package %smanifest", outType))
    out
}

################################################################################

### Convert the rgSet into the outType array
# Example usage:
# rgSetOutv1 <- combineArrays(rgSet1 = rgSetv1,
#                             rgSet2 = rgSetv2)
# input order does not matter here
.convertArray <- function(rgSet,
                         outType = c("IlluminaHumanMethylation450k",
                                     "IlluminaHumanMethylationEPIC",
                                     "IlluminaHumanMethylationEPICv2"),
                         verbose = verbose) {

    outType <- match.arg(outType)

    array <- BiocGenerics::annotation(rgSet)[["array"]]
    if(array == outType) stop("'rgSet' already in the 'outType' array type.")
    manifestIn <- .getManifest(array)
    manifestOut <- .getManifest(outType)
    annoOut <- if(outType == "IlluminaHumanMethylation450k"){
        "ilmn12.hg19"
    } else if(outType == "IlluminaHumanMethylationEPIC") {
        "ilm10b4.hg19"
    } else if(outType == "IlluminaHumanMethylationEPICv2") {
        "20a1.hg38"
    }

    keepAddresses <- list(
        I = NULL,
        II = NULL,
        SnpI = NULL,
        SnpII = NULL,
        Control = NULL)

    ### Probes of Type I
    probesIn <- getProbeInfo(manifestIn, type = "I")
    probesOut <- getProbeInfo(manifestOut, type = "I")
    if(array == "IlluminaHumanMethylationEPICv2") {
        # Remove duplicated probes
        probesIn$Name <- gsub("_.*", "", probesIn$Name)
        dupProbes <- table(probesIn$Name) %>%
            .[.>1] %>%
            names()
        probesIn <- probesIn[!probesIn$Name %in% dupProbes, ]
    }

    # Match the probes in the input and converted output
    commonNames <- intersect(probesIn$Name, probesOut$Name)
    probesIn <- probesIn[match(commonNames, probesIn$Name), ]
    probesOut <- probesOut[match(commonNames, probesOut$Name), ]

    # Remove probes with changed probe chemistry
    if(!all(probesIn$Color == probesOut$Color)){
        matchColor <- which(probesIn$Color == probesOut$Color)
        probesIn <- probesIn[matchColor, ]
        probesOut <- probesOut[matchColor, ]
    }

    # Sanity checks
    stopifnot(all(probesIn$Color == probesOut$Color))
    stopifnot(all(probesIn$ProbeSeqA == probesOut$ProbeSeqA))
    stopifnot(all(probesIn$ProbeSeqB == probesOut$ProbeSeqB))

    # Translating rgSetIn addresses to rgSetOut addresses
    translate <- c(probesOut$AddressA, probesOut$AddressB)
    names(translate) <- c(probesIn$AddressA, probesIn$AddressB)
    keepAddresses$I <- translate


    ### Probes of Type II
    probesIn <- getProbeInfo(manifestIn, type = "II")
    probesOut <- getProbeInfo(manifestOut, type = "II")
    if(array == "IlluminaHumanMethylationEPICv2") {
        # Remove duplicated probes
        probesIn$Name <- gsub("_.*", "", probesIn$Name)
        dupProbes <- table(probesIn$Name) %>%
            .[.>1] %>%
            names()
        probesIn <- probesIn[!probesIn$Name %in% dupProbes, ]
    }

    # Match the probes in the input and converted output
    commonNames <- intersect(probesIn$Name, probesOut$Name)
    probesIn <- probesIn[match(commonNames, probesIn$Name), ]
    probesOut <- probesOut[match(commonNames, probesOut$Name), ]

    # Remove probes with changed probe chemistry
    if(!all(probesIn$ProbeSeqA == probesOut$ProbeSeqA)){
        matchSeq <- which(probesIn$ProbeSeqA == probesOut$ProbeSeqA)
        probesIn <- probesIn[matchSeq, ]
        probesOut <- probesOut[matchSeq, ]
    }

    # Sanity checks
    stopifnot(all(probesIn$ProbeSeqA == probesOut$ProbeSeqA))

    # Translating rgSetIn addresses to rgSetOut addresses
    translate <- probesOut$AddressA
    names(translate) <- probesIn$AddressA
    keepAddresses$II <- translate


    ### Probes of Type SnpI
    probesIn <- getProbeInfo(manifestIn, type = "SnpI")
    probesOut <- getProbeInfo(manifestOut, type = "SnpI")
    if(array == "IlluminaHumanMethylationEPICv2") {
        # Remove duplicated probes
        probesIn$Name <- gsub("_.*", "", probesIn$Name)
        dupProbes <- table(probesIn$Name) %>%
            .[.>1] %>%
            names()
        probesIn <- probesIn[!probesIn$Name %in% dupProbes, ]
    }

    # Match the probes in the input and converted output
    commonNames <- intersect(probesIn$Name, probesOut$Name)
    probesIn <- probesIn[match(commonNames, probesIn$Name), ]
    probesOut <- probesOut[match(commonNames, probesOut$Name), ]

    # Sanity checks - sometimes the ProbeSeq A vs B is swapped - is that ok??
    stopifnot(all((probesIn$ProbeSeqA == probesOut$ProbeSeqB) |
                      (probesIn$ProbeSeqA == probesOut$ProbeSeqA)))
    stopifnot(all((probesIn$ProbeSeqB == probesOut$ProbeSeqB) |
                      (probesIn$ProbeSeqB == probesOut$ProbeSeqA)))

    # Translating rgSet2 addresses to rgSet1 addresses
    translate <- c(probesOut$AddressA, probesOut$AddressB)
    names(translate) <- c(probesIn$AddressA, probesIn$AddressB)
    keepAddresses$SnpI <- translate

    ### Probes of Type SnpII
    probesIn <- getProbeInfo(manifestIn, type = "SnpII")
    probesOut <- getProbeInfo(manifestOut, type = "SnpII")
    if(array == "IlluminaHumanMethylationEPICv2") {
        # Remove duplicated probes
        probesIn$Name <- gsub("_.*", "", probesIn$Name)
        dupProbes <- table(probesIn$Name) %>%
            .[.>1] %>%
            names()
        probesIn <- probesIn[!probesIn$Name %in% dupProbes, ]
    }

    # Match the probes in the input and converted output
    commonNames <- intersect(probesIn$Name, probesOut$Name)
    probesIn <- probesIn[match(commonNames, probesIn$Name),]
    probesOut <- probesOut[match(commonNames, probesOut$Name),]
    stopifnot(all(probesIn$ProbeSeqA == probesOut$ProbeSeqA))

    # Translating rgSet2 addresses to rgSet1 addresses
    translate <- probesOut$AddressA
    names(translate) <- probesIn$AddressA
    keepAddresses$SnpII <- translate

    ### Probes of Type Control
    probesIn <- getProbeInfo(manifestIn, type = "Control")
    probesOut <- getProbeInfo(manifestOut, type = "Control")
    commonAddress <- intersect(probesIn$Address, probesOut$Address)
    probesIn <- probesIn[match(commonAddress, probesIn$Address),]
    probesOut <- probesOut[match(commonAddress, probesOut$Address),]
    translate <- probesOut$Address
    names(translate) <- probesIn$Address
    keepAddresses$Control <- translate

    # Update rgSet

    keep <- unlist(unname(keepAddresses))
    rgSet <- rgSet[rownames(rgSet) %in% names(keep), ]
    rownames(rgSet) <- keep[rownames(rgSet)]
    annotation(rgSet) <- c(array = outType, annotation = annoOut)
    rgSet
}

################################################################################

### Combine multiple array types
# Example usage:
# rgSetv1 <- .convertArray(rgSetv2,
#                         outType = "IlluminaHumanMethylationEPIC",
#                         verbose = T)
.combineArrays <- function(rgSet1, rgSet2,
                           outType = NULL){

    # Check to see if we have sensible inputs
    stopifnot(is(rgSet1, "RGChannelSet") & is(rgSet2, "RGChannelSet"))

    # Obtain annotations - if output type is not specified,
    anno1 <- BiocGenerics::annotation(rgSet1)[["array"]]
    anno2 <- BiocGenerics::annotation(rgSet2)[["array"]]
    options <- c("IlluminaHumanMethylation450k",
                 "IlluminaHumanMethylationEPIC",
                 "IlluminaHumanMethylationEPICv2")
    if(is.null(outType)) {
        outType <- options[options %in% c(anno1, anno2)][1]
    } else {
        outType <- match.arg(outType, options)
    }
    message(paste0("combining data as ", outType, "\n"))

    if(anno1 == anno2) {
        features <- intersect(rownames(rgSet1),
                              rownames(rgSet2))
        rgSet <- BiocGenerics::combine(rgSet1[features, ],
                                       rgSet2[features, ])
    } else if(!anno1 %in% options | !anno2 %in% options) {
        stop("These array types cannot be combined.")
    } else {
        if(anno1 != outType) {
            rgSet1 <- .convertArray(rgSet1, outType = outType, verbose = verbose)
        }
        if(anno2 != outType) {
            rgSet2 <- .convertArray(rgSet2, outType = outType, verbose = verbose)
        }
        features <- intersect(rownames(rgSet1),
                              rownames(rgSet2))
        rgSet <- BiocGenerics::combine(rgSet1[features, ],
                                       rgSet2[features, ])
        rgSet$ArrayTypes <- rep(
            x = c(anno1, anno2),
            times = c(ncol(rgSet1), ncol(rgSet2)))
        rgSet
    }
}

################################################################################

.pickCompProbes2 <- function(betas, meta, nP, ct, ps = c("any", "both"),
                            trainingProbes = NULL, p.val = 1e-8,
                            min.delta.beta = 0, plot) {

    df <- as.matrix(betas)
    pd <- as.data.frame(meta)
    rownames(pd) <- pd$sampleNames

    ## only keep the cell type you're estimating
    keep <- which(pd$cellType %in% ct)
    pd <- pd[keep, ]
    df <- df[, keep]
    pd$cellType <- factor(pd$cellType, levels = ct) # make cell type a factor

    ## Find the probes with significant different methylation
    ## status across cell types
    tIndexes <- .splitit(pd$cellType)
    tIndexes <- tIndexes[lapply(tIndexes, length) > 0]
    if (!is.null(trainingProbes)) {
        trainingProbes <- trainingProbes[trainingProbes %in% rownames(df)]
        if (sum(trainingProbes %in% rownames(df)) / length(trainingProbes) < 0.9)
            message("Less than 90% of the training probes (likely IDOL) is
                    present in your dataset. This might impact prediction
                    accuracy.")
    } else {
        ## Identify which probe contribute to significantly different
        ## methylation status with t-test
        tstatList <- lapply(tIndexes, function(ct) {
            ctIndex <- rep(0, ncol(df))
            ctIndex[ct] <- 1
            return(genefilter::rowttests(df, factor(ctIndex)))
            # For each cell type, compare the mean DNAm level in this
            # cell type vs all others
        })
        ## Select N (default = 100) probes for each given cell type that
        ## can best distinguish cell types
        if (ps == "both"){
            probeList <- lapply(tstatList, function(ct) {
                probes <- ct[ct$p.value < p.val & abs(ct$dm) > min.delta.beta, ]
                pUp <- probes[order(probes$dm, decreasing = TRUE), ]
                pDown <- probes[order(probes$dm, decreasing = FALSE), ]
                return(c(rownames(pUp)[1:(nP/2)], rownames(pDown)[1:(nP/2)]))
            })
        } else if (ps == "any"){
            probeList <- lapply(tstatList, function(ct) {
                probes <- ct[ct$p.value < p.val & abs(ct$dm) > min.delta.beta, ]
                probes <- probes[order(abs(probes$dm), decreasing = TRUE), ]
                return(rownames(probes)[1:nP])
            })
        }
        trainingProbes <- unlist(probeList)
    }

    if (plot){
        pltct <- c(Bcell_cord = "#9C9EDEFF", Bnv = "#7375B5FF",
                   Bmem = "#4A5584FF", Bcell = "#7375B5FF",
                   CD4T_cord = "#CEDB9CFF", CD4nv = "#B5CF6BFF",
                   CD4mem = "#637939FF", Treg = "#8CA252FF", CD4T = "#B5CF6BFF",
                   CD8T_cord = "#E7CB94FF", CD8nv = "#E7BA52FF",
                   CD8mem = "#8C6D31FF", CD8T = "#E7BA52FF",
                   NK_cord = "#7BBCB0FF", NK = "#3A7C89FF",
                   Mono_cord = "#F3CBD3FF", Mono = "#707070",
                   Gran_cord = "#D39C83FF", Gran = "#A65461FF",
                   Neu = "#A65461FF", Bas = "#7B4173FF", Eos = "#A55194FF",
                   nRBC = "#843C39FF", PBMC = "#A5AA99", WBC = "#252525FF")
        pltct <- pltct[names(pltct) %in% ct]
        anncolors <- list(cellType = pltct)

        pheatmap::pheatmap(df[rownames(df) %in% trainingProbes, ],
                           show_colnames = F,
                           show_rownames = F,
                           annotation_col = pd[, "cellType", drop = F],
                           annotation_colors = anncolors)
    }


    ## Call a linear model with the selected probes and calculate the weights
    ## of each cell type
    df <- df[unique(trainingProbes), ]
    # Subset down the reference data to the selected probes -
    # remove duplicated CpG since they generate the same coefficients
    mod <- stats::model.matrix(~ pd$cellType - 1) %>% as.data.frame()
    colnames(mod) <- levels(pd$cellType)
    form <- stats::as.formula(sprintf("bt ~ %s - 1", paste(ct, collapse = " + ")))
    coef <- sapply(1:nrow(df), function(probe){
        mod$bt <- df[probe, ]
        fit <- stats::lm(form, data = mod[stats::complete.cases(mod), ])
        # Remove samples with missing methylation values
        fitCoef <- fit$coef
        return(fitCoef)
    }) %>% t() %>% as.data.frame()
    rownames(coef) <- rownames(df)

    return(coef)
}

################################################################################

.cp <- function(samp.n, coef, conditions = NULL) {

    requireNamespace("quadprog")
    nCt <- ncol(coef)
    nSamp <- ncol(samp.n)

    if(!is.null(conditions)) {
        message("Running on constraint")
        amat <- cbind(rep(-1, nCt), diag(nCt))
        b0 <- c(-1, rep(0, nCt))
    } else {
        message("No constraint")
        amat <- diag(nCt)
        b0 <- rep(0, nCt)
    }

    ctEst <- sapply(1:nSamp, function(samp){
        cpgs <- which(!is.na(samp.n[, samp]))
        dmat <- t(coef[cpgs, ]) %*% coef[cpgs, ]
        out <- quadprog::solve.QP(dmat,
                                  t(coef[cpgs, ]) %*% samp.n[cpgs, samp],
                                  amat,
                                  b0)$sol
        return(out)
    }) %>% t()
    rownames(ctEst) <- colnames(samp.n)
    colnames(ctEst) <- colnames(coef)

    return(ctEst)
}

################################################################################

.rpc <- function(samp.n, coef, conditions = 50) {

    ctEst <- apply(samp.n, 2, function(samp){
        fit <- MASS::rlm(x = coef, y = samp, maxit = conditions)
        out <- summary(fit)$coef[, 1]
        out[out < 0] <- 0
        out <- out/sum(out)
        # length(out)
        return(out)
    }) %>% t()
    return(ctEst)
}

################################################################################

.svr <- function(samp.n, coef, conditions = c(0.25, 0.5, 0.75)) {

    requireNamespace("e1071")
    ctEst <- lapply(conditions, function(nu){
        est <- apply(samp.n, 2, function(samp){
            fit <- e1071::svm(x = coef, y = samp,
                              scale = TRUE,
                              type = "nu-regression",
                       kernel = "linear", nu = nu)
            out <- t(fit$coefs) %*% fit$SV
            out[out < 0] <- 0
            out <- out/sum(out)
        }) %>% t()
        colnames(est) <- colnames(coef)
        rmse <- sqrt(colMeans((samp.n - coef %*% t(est))^2)) # calculate rmse
        return(list(est = est, rmse = rmse))
    })
    # select best nu
    rmse <- sapply(1:length(conditions), function(nv) ctEst[[nv]][["rmse"]])
    nuIndex <- apply(rmse, 1, which.min)
    ctEstF <- sapply(1:ncol(samp.n), function(samp){
        # select estimate for each sample based on best nu
        out <- ctEst[[nuIndex[samp]]]$est[samp, ]
        return(out)
    }) %>% t()
    rownames(ctEstF) <- colnames(samp.n)
    return(ctEstF)
}

################################################################################

.getErrorPerSample <- function(applyIndex, predictedIN,
                              coefDataIN, betasBulkIN) {
    trueBulk <- matrix(ncol = 1, nrow = nrow(coefDataIN),
                       data = 0)
    for (i in seq_len(ncol(coefDataIN))) {
        trueBulk[, 1] <- trueBulk[, 1] + coefDataIN[, i] *
            predictedIN[applyIndex, i]
    }
    betasBulkIN <- t(apply(betasBulkIN, 1, function(x) {
        x[is.na(x)] <- 0
        return(x)
    }))
    error <- .rmse(trueBulk, betasBulkIN[, applyIndex])
    return(error)
}

################################################################################

.hc2split <- function(x){
    A <- x$merge # (n-1,n) matrix
    n <- nrow(A) + 1
    B <- list()

    for (i in 1:(n-1)) {
        ai <- A[i, 1]

        if(ai < 0)
            B[[i]] <- -ai
        else
            B[[i]] <- B[[ai]]

        ai <- A[i, 2]

        if(ai < 0)
            B[[i]] <- sort(c(B[[i]], -ai))
        else
            B[[i]] <- sort(c(B[[i]], B[[ai]]))
    }

    CC <- matrix(rep(0, n*(n - 1)),
                 nrow = (n - 1),
                 ncol = n)

    for (i in 1:(n - 1)) {
        bi <- B[[i]]
        m <- length(bi)
        for(j in 1:m)
            CC[i,bi[j]] <- 1
    }

    split <- list(pattern = apply(CC, 1, paste, collapse = ""), member = B)

    return(split)
}

