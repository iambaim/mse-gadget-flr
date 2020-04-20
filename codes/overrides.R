# hcr_nafo is a function that produces a HCR that can have one or two stages. It has a Blim1 value (the one corresponding to precautionary approach Blim concept), below which F is zero. When SSB is between Blim1 and Btrigger1, F will be between 0 and Ftarget1; when SSB is bewteeen Btrigger1 and Blim2 (not related with the precautionary approach) F will be Ftarget1; when SSB is between Blim2 and Btrigger2, F will be between Ftarget1 and Ftarget2; and finally, when SSB is above Btrigger2, F will be equal to ftarget2. This would be a two stages HCR. If we want it only one stage (tipycal NAFO HCR, then Blim2 has to be defined as NULL)
nafo.hcr <- function (stk, fmin, ftarget1, blim1, btrigger1, ftarget2, blim2, btrigger2, ssb_lag, genArgs, tracking)
{
   ay <- genArgs$ay

   if(is.null(blim2)) {
      blim2=blim1
      btrigger2=btrigger1
      ftarget2=ftarget1
   }

   ssb <- ssb(stk)[, ac(ay - ssb_lag)]
   fout <- FLQuant(fmin, dimnames = list(iter = dimnames(ssb)$iter))
   fout[ssb >= btrigger2] <- ftarget2
   inbetween3 <- (ssb < btrigger2) & (ssb >= blim2)
   gradient3 <- (ftarget2 - ftarget1)/(btrigger2 - blim2)
   fout[inbetween3] <- (ssb[inbetween3] - blim2) * gradient3 + ftarget1
   inbetween2 <- (ssb < blim2) & (ssb >= btrigger1)
   fout[inbetween2] <- ftarget1
   inbetween1 <- (ssb < btrigger1) & (ssb >= blim1)
   gradient1 <- (ftarget1 - fmin)/(btrigger1 - blim1)
   fout[inbetween1] <- (ssb[inbetween1] - blim1) * gradient1 + fmin
   ctrl <- getCtrl(c(fout), "f", ay + 1, dim(fout)[6])
   list(ctrl = ctrl, tracking = tracking)
}


# New Stock recruitment function (true plus noise)
applyNoise.oem <- function (stk, idx, genArgs, stockName, ...)
{
    args <- list(...)

    # Get prev year
    ay <- genArgs$ay
    year <- ac(ay-1)    

    # Get noise parameters
    noiseParams.catch <- eval(parse(text=paste0(stockName, ".residual.params.catch")))
    noiseParams.index <- eval(parse(text=paste0(stockName, ".residual.params.index")))
    noiseParams.stock <- eval(parse(text=paste0(stockName, ".residual.params.stock")))

    saParam <- eval(parse(text=paste0(stockName, ".assessment")))

    print(paste("Applying noise before", saParam))

    # Multiply catch and catch.n
    #print("Before")
    #print(catch.n(idx[[1]])[,year])
    #print(index(idx[[1]])[,year])
    #print(fbar(stk)[,year])
    #print(ssb(stk)[,year])

    if(saParam == "SCAA" || saParam == "SAM") {
	if(!is.null(noiseParams.catch)) {
		noise.catch <- apply(noiseParams.catch, 1, function(x) exp(rnorm(1, x[["mean"]], x[["sd"]])))

		catch.n(stk)[,year] <- catch.n(stk)[,year] * noise.catch

		# Recalculate catch
		catch(stk)[,year] <- sum(catch.n(stk)[,year] * catch.wt(stk)[,year])

		# Recalculate harvest
		harvest(stk)[,year] <- -log((stock.n(stk)[,year] - catch.n(stk)[,year])/stock.n(stk)[,year])
	}
	if(!is.null(noiseParams.index)) {
		noise.index <- apply(noiseParams.index, 1, function(x) exp(rnorm(1, x[["mean"]], x[["sd"]])))

		multiplyIndex <- function(singleIdx, nsOEM) {
		    catch.n(singleIdx)[,year] <- catch.n(singleIdx)[,year] * nsOEM
		    index(singleIdx)[,year] <- index(singleIdx)[,year] * nsOEM
		    return(singleIdx)
		}

		idx <- lapply(idx, multiplyIndex, noise.index)
	}
    }else if(saParam == "truePlusNoise") {
	if(!is.null(noiseParams.stock)) {
		noise.stock <- apply(noiseParams.stock, 1, function(x) exp(rnorm(1, x[["mean"]], x[["sd"]])))

		stock.n(stk)[,year] <- stock.n(stk)[,year] * noise.stock

		# Recalculate stock
		stock(stk)[,year] <- sum(stock.n(stk)[,year] * stock.wt(stk)[,year])
	}
    }

    #print("After")
    #print(catch.n(idx[[1]])[,year])
    #print(index(idx[[1]])[,year])
    #print(fbar(stk)[,year])
    #print(ssb(stk)[,year])

    return (list(stk = stk, idx = idx))
}

# New Stock recruitment function (true plus noise)
truePlusNoise.sa <- function (stk, idx, ...) 
{
    args <- list(...)
    tracking <- args$tracking 
    tracking["conv.est", ac(range(stk)["maxyear"] + 1)] <- 0
    list(stk = stk, tracking = tracking)
}


tac.is.fixed <- function (stk, ctrl, genArgs, delta_tac_max = NA, delta_tac_min = NA,
    tracking) 
{
    iy <- genArgs$iy
    ay <- genArgs$ay
    it <- genArgs$it
    nsqy <- genArgs$nsqy
    if (ay == iy) 
        refCatch <- tracking["C.om", ac(ay - 1)]
    else refCatch <- tracking["metric.is", ac(ay - 1)]
    yrs <- as.numeric(dimnames(stock.n(stk))$year)
    last_data_yr <- yrs[length(yrs)]
    sqy <- (last_data_yr - nsqy + 1):last_data_yr
    fsq0 <- yearMeans(fbar(stk)[, ac(sqy)])
    ninter_yrs <- ay - last_data_yr
    ctrl <- getCtrl(c(rep(fsq0, times = ninter_yrs), ctrl@trgtArray[, 
        "val", ]), "f", (last_data_yr + 1):(ay + 1), dim(stock.n(stk))[6])
    nproj_yrs <- (ay + 1) - last_data_yr
    stkTmp <- stf(stk, nproj_yrs, wts.nyears = nsqy)
    gmean_rec <- c(exp(yearMeans(log(rec(stk)[, ac(sqy)]))))
    stkTmp <- fwd(stkTmp, ctrl = ctrl, sr = list(model = "mean", 
        params = FLPar(gmean_rec, iter = it)))
    TAC <- catch(stkTmp)[, ac(ay + 1)]
    upper_limit <- refCatch * delta_tac_max
    lower_limit <- refCatch * delta_tac_min
    TAC <- pmin(c(upper_limit), c(TAC), na.rm = TRUE)
    TAC <- pmax(c(lower_limit), c(TAC), na.rm = TRUE)
    ctrl <- getCtrl(c(TAC), "catch", ay + 1, it)
    list(ctrl = ctrl, tracking = tracking)
}

sam.sa <- function(stk, idx, ...)
{
    args <- list(...)
    stk.ctrl <- FLSAM.control(stk,idx)
    fit <- FLSAM(stk, idx, stk.ctrl, return.fit = TRUE)
    stk.sam <- SAM2FLR(fit, stk.ctrl)

    stk     <- stk + stk.sam

    tracking <- args$tracking

    # TODO: Check below
    if(attr(fit$opt$objective, "logarithm"))
        tracking["conv.est", ac(range(stk)["maxyear"] + 1)] <- -log(as.numeric(fit$opt$objective))
    else
        tracking["conv.est", ac(range(stk)["maxyear"] + 1)] <- as.numeric(fit$opt$objective)

    list(stk = stk, tracking = tracking)
}

