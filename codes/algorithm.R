############################################################################################################################
## Developing the FLa4a-Gadget MSE framework (IMR-REDUS project: 2018-2020) - with the haddock Gadget model
## started: mar 11, 2019; updated mar 12, 2019

#############################################################################################################################
## The a4a-gadget framework is designed to run MSE with a single or multispcies gadget model as an operating model
## 

## Implemented features
## 1. Multiple recruitment number can be generated from a matrix
## 2. Error (in what???) for the latest year is genarated using the function (of what???) 
## 3. MSEs can be run with parallel computing
## 4. MSEs can be run for 2 options using Switches:
##                1) Truth plus noise (perfect observation) - noise is applied to stock only.
##                2) Stock Assessment (w/ observation error) - noise is applied to both catch and index
## 5. SAM as an assessment model


## List of tasks to be done
## 1.
## 2.
## 3. 
##
##


##==============================================================================
## libraries and auxiliary functions
##==============================================================================
## Install latest gadgetr
#devtools::install_github("REDUS-IMR/gadget", ref="gadgetr")
## Install latest MSE from a4a

## load required packages
library(mse)
library(dplyr)
library(FLa4a)
library(FLash)
library(FLAssess)
library(ggplotFL)
library(FLBRP)
library(FLCore)
#library(MASS)
library(FLSAM)
library(filelock)


## set the global parameter values for MSE
frst.year <- 1988
prj.year <- 2017 #2017
fnl.year <- 2037 #2070
#
#
#


## a function to run MSE with 
## input parameters
##
##
runOneTimeline <- function(iterSim, saveRaw) {
	
	## Set the seed for random number generators
	set.seed(123)

	## Set directories
	codeDir <- paste0(homeDir,"/codes") # scripts the main MSE and auxiliary functions
	paramFileDir <- paste0(homeDir,"/paramfiles") # gadget model input parameters - REPLACE THIS FOLDER WHEN USING A DIFFERENT GADGET MODEL

	## load helper functions
	while("funcs:ExtraFunctions" %in% search()) detach("funcs:ExtraFunctions")
	sys.source(paste0(codeDir,"/funs.R"), attach(NULL, name = "funcs:ExtraFunctions"))

	## Load gadgetr locally
	library(gadgetr)

	## Load helper functions
	source(paste0(codeDir,"/gadget-fls.R"), local = T)
	source(paste0(codeDir,"/overrides.R"), local = T)

	##==============================================================================
	## Load the stock and gadget model parameter files
	##==============================================================================
	## set a directory with gadget model files
	#setwd("../models/LF_304")
	#setwd("../models/SC05_98_forecast")
	setwd(paste0(homeDir, "/models/", modelName)) # REPLACE THIS FOLDER WHEN USING A DIFFERENT GADGET MODEL

	## Load gadget model parameters
	#gadget(c("-s","-main","main.gadcap","-i","paramsin.gadcap"))
	gadget(c("-s","-main","main.gadcap","-i","paramsin.forecast"))

	## Initialize simulation
	initSim() # SOME DESCRIPTION HERE

	## list stocks from the Gadget model
	stockList <- c("cod", "red", "shrimp") # a lits of stocks in the operating model

	## specify stock-specific info from the Gadget model
	## STOCK 1: cod
	## STOCK 2: redfish
	## STOCK 3: shrimp
	
	## STOCK 1) cod ##############################################################################
	stk1.fleets <- c("cod.trawl_I", "cod.trawl_II", "cod.gil", "cod.long", "cod.trawl_forecast")
	stk1.stocks <- c("cod.hatching", "cod.imm", "cod.mat.small", "cod.mat.large")
	stk1.stocks.mature <- c("cod.mat.small", "cod.mat.large")
	stk1.surveys <- c("cod.surveyEU", "cod.surveyEU_forecast")
	stk1.forecasts <- c("cod.trawl_forecast")
	stk1.forecasts.tac.proportion <- c(0.232, 0.351, 0.298, 0.119) # SOME DESRIPTIUON HERE - where did these number come from???
	
	## specify reference points
	stk1.hcr.params <- list(method=nafo.hcr, args=list(ssb_lag=1, 
	                                                  fmin=0.0000001, 
	                                                  blim1=17906000, 
	                                                  btrigger1=25943000, 
	                                                  ftarget1=fComb[combIndex,"cod"], 
	                                                  blim2=45000000, 
	                                                  btrigger2=55000000, 
	                                                  ftarget2=0.55))
	
	## m2=NULL: we calculate m2 from gadget result
	## m2=0 means we use only residual mortality (m1)
	stk1.params <- list(minage=1, 
	                   maxage=12, 
	                   minfbar=3, 
	                   maxfbar=7, 
	                   startf=0.56, 
	                   endf=0.65, 
	                   m1=c(0.35), m2=NULL)
	
	## Recruitment parameters, if read csv (data frame) will apply the values accordingly, if a constant value, 
	## will apply the value as mux, if NULL leaving the recruitment params as it is
	stk1.recruit.params <- read.csv(paste0(paramFileDir, "/muxfactors_cod.csv"))
	stk1.recruit.params <- stk1.recruit.params[, c(1, iterSim + 1)]
	#cod.recruit.params <- 10.75871423
	
	## there are 2 options for MP (1. truePlusNoise and 2. SCAA)
	## If truePlusNoise is chosen the noise using the residual params will be applied to stock only (perfect observation)
	## If sca.sa is chosen, the noise will be applied to both catch and index
	stk1.residual.params.catch <- NULL #read.csv(paste0(paramFileDir, "/cod_resid_pars_catch.csv"))
	stk1.residual.params.index <- NULL #read.csv(paste0(paramFileDir, "/cod_resid_pars_index.csv"))
	stk1.residual.params.stock <- NULL #read.csv(paste0(paramFileDir, "/cod_resid_pars_stock.csv"))
	stk1.assessment <- "truePlusNoise" #truePlusNoise or SCAA
	## If you don't want to apply error:
	#stk1.residual.params <- NULL
	stk1.noteating.forecast <- TRUE

	# ## STOCK 2) redfish ###############################################################################
	# red.fleets <- c("red.trawl_I", "red.trawl_II", "shrimp.red.trawl", "red.trawl_forecast")
	# red.stocks <- c("red.hatching", "red.fem.imm", "red.fem.matu", "red.male.imm", "red.male.matu")
	# red.stocks.mature <- c("red.fem.matu")
	# red.surveys <- c("red.surveyEU", "red.surveyEU_forecast")
	# red.forecasts <- c("red.trawl_forecast")
	# red.forecasts.tac.proportion <- c(0.232, 0.351, 0.298, 0.119)
	# 
	# ## specify reference points
	# red.hcr.params <- list(method=nafo.hcr, 
	#                        args=list(ssb_lag=1, 
	#                                  fmin=0.0000001, 
	#                                  blim1=22027000, 
	#                                  btrigger1=35361000, 
	#                                  ftarget1=fComb[combIndex,"red"], 
	#                                  blim2=NULL, 
	#                                  btrigger2=NULL, 
	#                                  ftarget2=NULL))
	# 
	# ## m2=NULL means we calculate m2 from gadget result
	# ## m2=0 means we use only residual mortality (m1)
	# red.params <- list(minage=1, 
	#                    maxage=25, 
	#                    minfbar=6, 
	#                    maxfbar=16, 
	#                    startf=0.56, 
	#                    endf=0.65, 
	#                    m1=c(0.145, 0.101, 0.072, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
	#                         0.05, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.12, 
	#                         0.15, 0.2, 0.25, 0.25, 0.25, 0.25), 
	#                    m2=NULL)
	# 
	# ## Recruitment parameters, if read csv (data frame) will apply the values accordingly, if a constant value, 
	# ## will apply the value as mux, if NULL leaving the recruitment params as it is	
	# red.recruit.params <- read.csv(paste0(paramFileDir, "/muxfactors_red.csv"))
	# red.recruit.params <- red.recruit.params[, c(1, iterSim+1)]
	# #red.recruit.params <- 33.23606361
	# 
	# ## 2 options for an assessment model (1. truePlusNoise and 2. SCAA)
	# ## If truePlusNoise is chosen the noise using the residual params will be applied to stock only
	# ## If sca.sa is chosen, the noise will be applied to both catch and index
	# red.residual.params.catch <- NULL #
	# red.residual.params.index <- NULL #
	# red.residual.params.stock <- NULL #
	# red.assessment <- "truePlusNoise" # truePlusNoise or SCAA
	# red.noteating.forecast <- TRUE
	# 
	# ## STOCK 3) shrimp ######################################################################################## 
	# shrimp.fleets <- c("shrimp.trawl", "shrimp.trawl_forecast")
	# shrimp.stocks <- c("shrimp.hatching", "shrimp.fem.multi", "shrimp.fem.primi", "shrimp.male")
	# shrimp.stocks.mature <- c("shrimp.fem.multi")
	# shrimp.surveys <- c("shrimp.surveyEU", "shrimp.surveyEU_forecast")
	# shrimp.forecasts <- c("shrimp.trawl_forecast")
	# shrimp.forecasts.tac.proportion <- c(0.232, 0.351, 0.298, 0.119)
	# 
	# ## refrence points - derived from NAFO HCR
	# shrimp.hcr.params <- list(method=nafo.hcr, 
	#                           args=list(ssb_lag=1, 
	#                                     fmin=0.0000001, 
	#                                     blim1=11864000, 
	#                                     btrigger1=31114000, 
	#                                     ftarget1=fComb[combIndex,"shrimp"], 
	#                                     blim2=NULL, 
	#                                     btrigger2=NULL, 
	#                                     ftarget2=NULL))
	# 
	# ## m2=NULL: we calculate m2 from gadget result
	# ## m2=0: we use only residual mortality (m1)
	# shrimp.params <- list(minage=1, 
	#                       maxage=7, 
	#                       minfbar=2, 
	#                       maxfbar=5, 
	#                       startf=0.56, 
	#                       endf=0.65, 
	#                       m1=c(0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), 
	#                       m2=NULL)
	# shrimp.recruit.params <- read.csv(paste0(paramFileDir, "/muxfactors_shrimp.csv"))
	# shrimp.recruit.params <- shrimp.recruit.params[, c(1, iterSim+1)]
	# #shrimp.recruit.params <- 2765.125525
	# 
	# ## 2 options for assessment models (1. truePlusNoise and 2. SCAA)
	# ## If truePlusNoise is chosen the noise using the residual params will be applied to stock only
	# ## If sca.sa is chosen, the noise will be applied to both catch and index
	# shrimp.residual.params.catch <- NULL #read.csv(paste0(paramFileDir, "/shrimp_resid_pars_catch.csv"))
	# shrimp.residual.params.index <- NULL #read.csv(paste0(paramFileDir, "/shrimp_resid_pars_index.csv"))
	# shrimp.residual.params.stock <- NULL #read.csv(paste0(paramFileDir, "/shrimp_resid_pars_stock.csv"))
	# shrimp.assessment <- "truePlusNoise" #truePlusNoise or SCAA
	# shrimp.noteating.forecast <- TRUE

	
	## set global parameters for the Gadget simulations
	firstYear <- frst.year
	projYear <- prj.year
	finalYear <- fnl.year
	
	## initialize gadget forecasting output
	gadgetOut <-list()

	## Run (WHAT?) until the first projection year
	gadcapOut <- runUntil(projYear-1)

	## Set the MSE loop parameters for each stock
	prepareStock  <- function(stockNameGl) {

		## specify a stock
		gadcap.ret <- gadcapOut[[stockNameGl]]
		stk <- gadcap.ret$stk
		idx <- FLIndices(a=gadcap.ret$idx)

		##==============================================================================
		## set MSE parameters
		##==============================================================================
    ## NEED A LOT MORE DESCRIPTION HERE
		it <- 1 # iterations
		fy <- finalYear #final year
		y0 <- range(stk)["minyear"] # initial data year
		dy <- range(stk)["maxyear"] # final data year
		iy <- projYear # initial year of projection (also 'intermediate' year)
		ny <- fy - iy + 1 # number of years to project from intial year
		nsqy <- 3 # number of years to compute status quo metrics
		vy <- ac(iy:fy) # vector of years to be projected

		## compute biological parameters of the stock - based on means of the last 5 years of observations
		stk <- stf(stk, fy-dy, nsqy, nsqy)

		##==============================================================================
		## Fleet behaviour
		##==============================================================================
		## NEED A LOT MORE DESCRIPTION HERE
		fb <- mseCtrl(method=hyperstability.fb, 
		              args=list(beta=0.8)) # WHAT IS BETA?

		##==============================================================================
		## OM object
		##==============================================================================
		## NEED A LOT MORE DESCRIPTION HERE
		om <- FLom(stock=stk)#, 
		        #fleetBehaviour=fb)
		#save(om, it, fy, y0, dy, iy, ny, nsqy, vy, fit, file="om.RData")

		###############################################################################
		## OEM settings
		###############################################################################
    ## SOME DESCRIPTION HERE OR DELETE
		
		
		##==============================================================================
		## prepare mse objects
		##==============================================================================
    ## NEED A LOT MORE DSCRIPTION HERE
		idx <- FLIndices(a=gadcap.ret$idx) # survey indices
		stk <- stock(om) # stock OM
		stk0 <- stk # initial stock object?

		##==============================================================================
		## compute the indices catchability from the a4a fit (without simulation)
		##==============================================================================

		## initialize indices
		idcs <- FLIndices()
		for (i in 1:length(idx)){
		  
			## assume surveys are conducted on Jan. 1st
			lst <- mcf(list(idx[[i]]@index, stock.n(stk0)))
			
			## log-transform catchability of the index
			idx.lq <- log(lst[[1]]/lst[[2]])
			
			## initialize stock objects
			idx.qmu <- idx.qsig <- stock.n(iter(stk,1))
			
			## catchability 
			## in this model catchability is constant
			idx.qmu[] <- yearMeans(idx.lq)
			idx.qsig[] <- sqrt(yearVars(idx.lq))
			idx.q <- FLQuant(NA, dimnames=dimnames(stock.n(stk)))
			
			## compute index catchability with noise generated from the lognormal distribution 
			## (mean and sd calculated above)
			idx.q <- rlnorm(it, idx.qmu, idx.qsig)
			#idx.q[,ac(y0:iy)] <- idx.q[,ac(y0:iy)]
			idx_temp <- idx.q * stock.n(stk)
			
			## generate initial index
			idx_temp <- FLIndex(index=idx_temp, 
			                    index.q=idx.q)
			range(idx_temp)[c("startf", "endf")] <- c(0, 0)
			idcs[[i]] <- idx_temp
		}
		
		names(idcs) <- names(idx)
		#idx <- FLIndices(a=idcs$a)

		##==============================================================================
		## Deviances for catch.n
		##==============================================================================
		## SOME DECRITPION HERE
		#catch.dev <- log(catch.n(stk))
		#catch.dev <- catch.dev-iterMeans(catch.dev)
		#Sig <- apply(catch.dev[,ac(y0:dy),1,1,,drop=TRUE], 3, function(x) cov(t(x)))
		#Sig <- apply(Sig, 1, mean)
		#Sig <- matrix(Sig, ncol=dim(catch.dev)[1])
		#catch.dev[,ac(vy)][] <- t(mvrnorm(it * length(vy), rep(0, nrow(Sig)), Sig))
		#catch.dev <- exp(catch.dev)

		##==============================================================================
		## OEM object
		##==============================================================================
    ## SOME DESCRITPOIN HERE
		idxDev <- lapply(idcs, index.q) # compute noise using the index catchbilaity computed above
		names(idxDev) <- "index.q"
		stkDev <- FLQuant()
		dev <- list(idx=idxDev, stk=stkDev)
		obs <- list(idx=idcs[1], stk=stk)
		#oem <- FLoem(method=sampling.oem, args=list(oe="index"), observations=obs, deviances=dev)
		oem <- FLoem()
		#save(oem, file="oem.RData")

		###############################################################################
		## Implementation error
		###############################################################################
    ## SOME DESCRITPOIN HERE
		#iem <- FLiem(method=noise.iem, args=list(fun="rlnorm", mean=0, sd=0.1, multiplicative=TRUE))
		iem <- FLiem()

		###############################################################################
		## Management procedure
		###############################################################################
    ## SOME DESCRIPTION HERE
		# specify general pars
		mpPars <- list(seed=1234, 
		               fy=fy, 
		               y0=y0, 
		               iy=iy, 
		               nsqy=nsqy, 
		               it=it)

		#==============================================================================
		## Scenarios
		#==============================================================================
		## Tell stocks to stop eating (if requested) - FOR PREDATION/FOOD CONSUMPTION IN GADGET?????
		if(eval(parse(text=paste0(stockNameGl, ".noteating.forecast")))){
		  stockCat <- eval(parse(text=paste0(stockNameGl, ".stocks")))
			tmp <- lapply(stockCat, stopEating)
			print("Stocks stop eating now")
			print(tmp)
			}

		## Get HCR parameters  NEED MORE SPECIFIC DESCRIPTION HERE
		hcrParams <- eval(parse(text=paste0(stockNameGl, ".hcr.params")))

		## Get stock assessment model parameter NEED MORE SPECIFIC DESCRIPTION HERE
		saParam <- eval(parse(text=paste0(stockNameGl, ".assessment")))
		if(saParam == "truePlusNoise") saMethod <- truePlusNoise.sa
		else if(saParam == "SCAA") saMethod <- sca.sa
		else if(saParam == "SAM") saMethod <- sam.sa

		## baseline model with TAC
		ctrl <- list(ctrl.hcr = mseCtrl(method=hcrParams[["method"]], 
		                                args=hcrParams[["args"]]),
			ctrl.is = mseCtrl(method=tac.is.fixed),
			ctrl.sa = mseCtrl(method=saMethod))

		## specify alternative scenarios
		scenarioName <- paste0(stockNameGl, ".", "iter", iterSim)

		return(list(opModel=om, 
		            indices=idx, 
		            obsModel=oem, 
		            impModel=NULL, 
		            ctrl.mp=ctrl, 
		            mpPars=mpPars, 
		            scenario=scenarioName, 
		            tracking=NULL))
	}

	## load Gadget helper functions
	source(paste0(codeDir,"/gadget-fwd.R"), local=T)
	source(paste0(codeDir,"/mp-methods-gadget.R"), local=T)
	inputPre <- lapply(stockList, prepareStock)
	names(inputPre) <- stockList
	res <- mp.gadget(inputPre)

	#return(list(mseResults=res,gadgetResults=gadgetOut))
	return(list(mseResults=res))
}

## specify MSE run parameters 
## Enable below to run directly from R shell
combIndex <- 1 # the index value for ???
iterIndex <- 1 # the number of iterations?

## set global Gadget model variables
homeDir <- paste0(getwd()) #"/../")
modelName <- "SC05_98_forecast_new_survey_forecast"
saveAllRawData <- FALSE

## Read effort combination
print(paste("I'm running with combination no.", combIndex, 
            "iteration", iterIndex))
fComb <- read.csv(paste0(homeDir, "/paramfiles/effort_combination.csv"))

## Run with combination and iterIndex
resultFinal <- runOneTimeline(iterIndex, saveAllRawData)

## set a directory for and save the output as a .rds file 
outFileName <- paste0(homeDir,"/results-combination", combIndex,".rds")

## Use lock to prevent race condition when combining results
lck <- lock(paste0(outFileName,".lock"))

## Check if there is any existing out put file in the folder
if(file.exists(outFileName)) {
	
  ## Load old results
	allResults <- readRDS(outFileName)
	} else {
	## initialize the output if there is no existing output file
	allResults <- list()
	}

## Combine all the outputs 
allResults[[iterIndex]] <- resultFinal

## Save the combined output info
write.table(fComb[combIndex,], 
            file=(paste0(outFileName,".info.txt")))

## Save all the outputs as a .rds file
saveRDS(allResults, file=outFileName)

## Unlock the file
unlock(lck)


