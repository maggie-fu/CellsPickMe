##############################################################################################################

splitit <- function(x) {
    split(seq(along = x), x)
}

##############################################################################################################


pickCompProbes2 <- function(betas, meta, nP, ct, ps = "any", trainingProbes = NULL, p.val = 1e-8, min.delta.beta = 0, plot) {

    df <- as.matrix(betas)
    pd <- as.data.frame(meta)
    rownames(pd) <- pd$sampleNames

    ## only keep the cell type you're estimating
    keep <- which(pd$cellType %in% ct)
    pd <- pd[keep, ]
    df <- df[, keep]
    pd$cellType <- factor(pd$cellType, levels = ct) # make cell type a factor

    ## Find the probes with significant different methylation status across cell types
    tIndexes <- splitit(pd$cellType)
    tIndexes <- tIndexes[lapply(tIndexes, length) > 0]
    if (!is.null(trainingProbes)) {
        trainingProbes <- trainingProbes[trainingProbes %in% rownames(df)]
        if (sum(trainingProbes %in% rownames(df)) / length(trainingProbes) < 0.9)
            message("Less than 90% of the training probes (likely IDOL) is present in your dataset. This might impact prediction accuracy.")
    } else {
        ## Identify which probe contribute to significantly different methylation status with t-test
        tstatList <- lapply(tIndexes, function(ct) {
            ctIndex <- rep(0, ncol(df))
            ctIndex[ct] <- 1
            return(genefilter::rowttests(df, factor(ctIndex))) # For each cell type, compare the mean DNAm level in this cell type vs all others
        })
        ## Select N (default = 100) probes for each given cell type that can best distinguish cell types
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
        pltct <- c(Bcell_cord = "#9C9EDEFF", Bnv = "#7375B5FF", Bmem = "#4A5584FF",
                   CD4T_cord = "#CEDB9CFF", CD4nv = "#B5CF6BFF", CD4mem = "#637939FF", Treg = "#8CA252FF",
                   CD8T_cord = "#E7CB94FF", CD8nv = "#E7BA52FF", CD8mem = "#8C6D31FF",
                   NK_cord = "#7BBCB0FF", NK = "#3A7C89FF", Mono_cord = "#F3CBD3FF", Mono = "#707070",
                   Gran_cord = "#D39C83FF", Neu = "#A65461FF", Bas = "#7B4173FF", Eos = "#A55194FF",
                   nRBC = "#843C39FF", PBMC = "#A5AA99", WBC = "#252525FF")
        pltct <- pltct[names(pltct) %in% ct]
        anncolors <- list(cellType = pltct)

        pheatmap::pheatmap(df[rownames(df) %in% trainingProbes, ],
                           show_colnames = F,
                           show_rownames = F,
                           annotation_col = pd[, "cellType", drop = F],
                           annotation_colors = anncolors)
    }


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
}

#############################################################################################################

CP <- function(samp.n, coef, conditions = NULL) {

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
        out <- quadprog::solve.QP(dmat, t(coef[cpgs, ]) %*% samp.n[cpgs, samp], amat, b0)$sol
        return(out)
    }) %>% t()
    rownames(ctEst) <- colnames(samp.n)
    colnames(ctEst) <- colnames(coef)

    return(ctEst)
}


#############################################################################################################

RPC <- function(samp.n, coef, conditions = 50) {

    ctEst <- apply(samp.n, 2, function(samp){
        fit <- MASS::rlm(x = coef, y = samp, maxit = conditions)
        out <- summary(fit)$coef[, 1]
        out[out < 0] <- 0
        out <- out/sum(out)
        length(out)
        # return(out)
    }) %>% t()
    return(ctEst)
}


#############################################################################################################


SVR <- function(samp.n, coef, conditions = c(0.25, 0.5, 0.75)) {

    requireNamespace("e1071")
    ctEst <- lapply(conditions, function(nu){
        est <- apply(samp.n, 2, function(samp){
            fit <- e1071::svm(x = coef, y = samp, scale = TRUE, type = "nu-regression",
                       kernel = "linear", nu = nu)
            out <- t(fit$coefs) %*% fit$SV
            out[out < 0] <- 0
            out <- out/sum(out)
        }) %>% t()
        colnames(est) <- colnames(coef)
        rmse <- sqrt(colMeans((samp.n - coef %*% t(est))^2)) # calculate rmse
        return(list(est = est, rmse = rmse))
    })
    rmse <- sapply(1:length(conditions), function(nv) return(ctEst[[nv]][["rmse"]])) # select best nu
    nuIndex <- apply(rmse, 1, which.min)
    ctEstF <- sapply(1:ncol(samp.n), function(samp){ # select estimate for each sample based on best nu
        out <- ctEst[[nuIndex[samp]]]$est[samp, ]
        return(out)
    }) %>% t()
    rownames(ctEstF) <- colnames(samp.n)
    return(ctEstF)
}


#############################################################################################################

getErrorPerSample <- function(applyIndex, predictedIN,
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
    error <- RMSE(trueBulk, betasBulkIN[, applyIndex])
    return(error)
}


#############################################################################################################

hc2split <- function(x){
    A <- x$merge # (n-1,n) matrix
    n <- nrow(A) + 1
    B <- list()

    for(i in 1:(n-1)){
        ai <- A[i,1]

        if(ai < 0)
            B[[i]] <- -ai
        else
            B[[i]] <- B[[ai]]

        ai <- A[i,2]

        if(ai < 0)
            B[[i]] <- sort(c(B[[i]],-ai))
        else
            B[[i]] <- sort(c(B[[i]],B[[ai]]))
    }

    CC <- matrix(rep(0,n*(n-1)),nrow=(n-1),ncol=n)

    for(i in 1:(n-1)){
        bi <- B[[i]]
        m <- length(bi)
        for(j in 1:m)
            CC[i,bi[j]] <- 1
    }

    split <- list(pattern=apply(CC,1,paste,collapse=""), member=B)

    return(split)
}

