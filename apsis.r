#!/usr/bin/R

## APSIS: Analysis of Power in Sequencing-and-Imputation Studies

library(data.table)

apsis_rds <- 'apsis_dat.rds'

apsis_dat <- readRDS(apsis_rds)

diseaseModel <- function(..., model_obj = NULL){
	if(is.null(model_obj)){
		model_obj <- list(...)
		if( all(c('model','mafs') %in% names(model_obj[[1]])) ){
			invisible(return(model_obj[[1]]))
		}
	}
	if( is.null(model_obj$freqs) ){
		names(model_obj) <- tolower(names(model_obj))
		
		kp_pars <- c('maf', 'rr', 'beta', 'prev', 'type')
		
		hwe <- function(x) c((1-x)^2, 2*x*(1-x), x^2)
		mom <- function(pmf, pow = 1) sum(((0:2)^pow)*pmf)
		nml <- function(x){
			if(any(x < 0) ) cat("\nWARNING: negative frequency coerced to 0\n")
			x*(x>=0)/sum(x*(x>=0))
		}
		
		model_text <- c()
		freqs <- list()
		
		# if( !all( names(model_obj) %in% kp_pars ) ){
			# cat("\nValid disease model parameters: 'maf', 'rr', 'beta', 'prev', 'type'\n\n")
			# stop(paste("Unknown parameters specified:",paste(setdiff(names(model_obj),kp_pars),collapse=',')))
		# }
		
		if( is.null(model_obj$type) ) model_obj$type <- ""
		
		if( is.null(model_obj$maf) ){
			stop('MAF not specified in disease model.')
		}else{
			model_obj$maf <- as.numeric(model_obj$maf)
			maf <- model_obj$maf
			model_text <- c(model_text, paste0('MAF=',round(maf*100,2),'%'))		
			if( !( maf > 0 & maf <= 0.5) ) stop("Please specify 0 < MAF <= 0.5 in disease model.")
			if( maf < 5e-04 ) cat("Warning: Estimates for very low MAF may be unreliable.")
			freqs$population <- nml(hwe(maf))
		}
		
		if( is.null(model_obj$prev) ){
			stop('Disease prevalence ("prev") not specified in disease model.')
		}else{
			model_obj$prev <- as.numeric(model_obj$prev)
			prev <- model_obj$prev
		}
		
		if( !is.null(model_obj$rr) ){
			model_obj$rr <- as.numeric(model_obj$rr)
			model_obj$beta <- NA
			rr <- model_obj$rr
			model_text <- c(model_text, paste0('rr=',rr))
			if( tolower(model_obj$type)=="recessive" ){
				freqs$cases <- c(NA,NA,rr*maf*maf/((rr-1)*maf*maf+1))
				freqs$cases[1] <- (1-maf)*(1-freqs$cases[3])/(1+maf)
				freqs$cases[2] <- 2*maf*(1-freqs$cases[3])/(1+maf)
			}else if( tolower(model_obj$type)=="dominant" ){
				freqs$cases <- c(((1-maf)^2)/(rr+(1-rr)*((1-maf)^2)), NA, NA)
				freqs$cases[3] <- maf*(1-freqs$cases[1])/(2 - maf)
				freqs$cases[2] <- 2*(1-maf)*freqs$cases[3]/maf
			}else if( model_obj$type %in% c("","additive","multiplicative") ){
				freqs$cases  <- hwe(rr*maf / (1 - maf + rr*maf))
			}else{
				stop(paste0('Unknown model type "',model_obj$type,'"; please specify "dominant", "recessive", or "additive" (default)'))
			}
			freqs$controls <- (freqs$population - prev*freqs$cases)/(1-prev)
		}else if( !is.null(model_obj$beta) ){
			model_obj$beta <- as.numeric(model_obj$beta)
			model_obj$rr <- NA
			model_text <- c(model_text, paste0('beta=',model_obj$beta))
			if( tolower(model_obj$type)=="recessive" ){
				betas <- c(0,0,1)*model_obj$beta
			}else if( tolower(model_obj$type)=="dominant" ){
				betas <- c(0,1,1)*model_obj$beta
			}else if( model_obj$type %in% c("","additive","multiplicative") ){
				betas <- c(0,1,2)*model_obj$beta
			}else{
				stop(paste0('Unknown model type "',model_obj$type,'"; please specify "dominant", "recessive", or "additive" (default)'))
			}
			E_g <- sum(freqs$population*betas)
			V_g <- sum(freqs$population*(betas^2)) - E_g^2
			thresh <- qnorm(1 - prev)
			EY <- 1 - pnorm((thresh - (betas - E_g) ) / sqrt(1 - V_g))
			EYn <- EY*prev/sum(EY*freqs$population)
			freqs$cases <- EYn*freqs$population/prev
			freqs$controls <- (freqs$population - prev*freqs$cases)/(1-prev)
		}else{
			stop('Please specify "rr" (relative risk) or "beta" (liability threshold) in disease model.')
		}
		
		model_text <- c(model_text, paste0('Prev=',round(prev*100,2),'%'))
		if(model_obj$type!="") model_text <- c(model_text, model_obj$type)
		
		freqs$cases <- nml(freqs$cases)
		freqs$controls <- nml(freqs$controls)
		
		invisible(list('model'=model_obj, 'freqs'=freqs, 'mafs' = lapply(freqs,mom), 'model_text' = paste(model_text, collapse = ',')))
	}else{
		model_obj
	}
}

diseaseModelList <- function(...){
	inp <- expand.grid(...)
	out <- list()
	for(i in 1:nrow(inp)){
		out[[i]] <- diseaseModel(model_obj=as.list(inp[i,]))
	}
	out$model_list <- TRUE
	invisible(out)
}

caseControlNCP <- function(PrG_Y1, PrG_Y0, w_seq = 0.5){
	mom <- function(pmf, pow = 1) sum(((0:2)^pow)*pmf)

	EG_Y1 <- mom(PrG_Y1); EG_Y0 <- mom(PrG_Y0)
	VG_Y1 <- mom(PrG_Y1,2) - EG_Y1^2
	VG_Y0 <- mom(PrG_Y0,2) - EG_Y0^2
	VG <- w_seq*VG_Y1 + (1 - w_seq)*VG_Y0 + ((EG_Y1 - EG_Y0)^2)*w_seq*(1 - w_seq)
	CGY <- w_seq*(1-w_seq)*(EG_Y1 - EG_Y0)
	
	CGY / sqrt(VG * w_seq*(1-w_seq))
}

getStudyMAF <- function(disease_model, cc_f){
	if(is.null(disease_model$freqs)) disease_model <- diseaseModel(disease_model)
	with(disease_model$mafs,
		cases*cc_f + controls*(1-cc_f)
	)
}

getNCP <- function(disease_model, cc_f){
	if(is.null(disease_model$freqs)) disease_model <- diseaseModel(disease_model)
	caseControlNCP( disease_model$freqs$cases, disease_model$freqs$controls, cc_f )
}

checkPopName <- function(population){
	if( population %in% names(apsis_dat$pop_labels) ){
		population
	}else{
		if( population %in% apsis_dat$pop_labels ){
			names(apsis_dat$pop_labels)[population==apsis_dat$pop_labels]
		}else{
			stop(paste("Please specify a valid population:\n",
				paste(apsis_dat$pop_labels,paste0('(',names(apsis_dat$pop_labels),')'),collapse=', ')
			))
		}
	}
}

checkArrayName <- function(array_name){
	if( array_name %in% names(apsis_dat$array_labels) ){
		array_name
	}else{
		if( array_name %in% apsis_dat$array_labels ){
			names(apsis_dat$array_labels)[array_name==apsis_dat$array_labels]
		}else{
			stop(paste("Please specify a valid genotyping array:\n",
				paste(apsis_dat$array_labels,paste0('(',names(apsis_dat$array_labels),')'),collapse=', ')
			))
		}
	}
}

apsisPower <- function(disease_model, n_imputed, n_sequenced, case_control_ratio = 1, population = "SR", snp_array = "OmniExpress", rsq_thresh = 0.3, alpha_level = 5e-8){
	
	if( !is.null(disease_model$model_list) ){
		if( disease_model$model_list ){
			disease_model$model_list <- NULL
			return(rbindlist(lapply(disease_model, apsisPower, n_sequenced=n_sequenced, n_imputed=n_imputed, case_control_ratio=case_control_ratio,population=population,snp_array=snp_array,rsq_thresh=rsq_thresh,alpha_level=alpha_level)))
		}
	}
	
	if( length(n_sequenced) > 1 ){
		return(rbindlist(lapply(n_sequenced, apsisPower, disease_model=disease_model, n_imputed=n_imputed, case_control_ratio=case_control_ratio,population=population,snp_array=snp_array,rsq_thresh=rsq_thresh,alpha_level=alpha_level)))
	}
	
	if( length(population) > 1 ){
		return(rbindlist(lapply(population, apsisPower, disease_model=disease_model, n_imputed=n_imputed, n_sequenced=n_sequenced, case_control_ratio=case_control_ratio,snp_array=snp_array,rsq_thresh=rsq_thresh,alpha_level=alpha_level)))
	}
	
	if( length(snp_array) > 1 ){
		return(rbindlist(lapply(snp_array, apsisPower, disease_model=disease_model, n_imputed=n_imputed,n_sequenced=n_sequenced, case_control_ratio=case_control_ratio,population=population,rsq_thresh=rsq_thresh,alpha_level=alpha_level)))
	}
	
	r <- case_control_ratio
	N <- n_sequenced
	M <- n_imputed
	
	POP <- checkPopName(population)
	POP_NAME <- apsis_dat$pop_labels[POP]
	ARR <- checkArrayName(snp_array)
	ARR_NAME <- apsis_dat$array_labels[ARR]
	ARR_COST <- apsis_dat$array_prices[ARR_NAME]
	SEQ_COST <- 1000
	disease_model <- diseaseModel(disease_model)
	
	eta <- getNCP(disease_model, r/(r + 1))

	N_total <- apsis_dat$pop_s_sizes[POP]
	
	N_max <- ifelse(POP=='LA', 1500, 2000)
	
	pr_pspec_var <- 'p.smooth'
	
	if( N > N_max ){
		stop(paste("Please specify n_sequenced <=",N_max,"for",POP_NAME))
	}
	
	c_a <- qnorm(0.5 * alpha_level, low = FALSE)
	
	maf_gwas <- getStudyMAF(disease_model, r/(r + 1))
	maf_pop <- disease_model$mafs$population
	
	max.m <- max(apsis_dat$RsqList[[POP]][[ARR]][[round(N/100)+1]]$m, na.rm = TRUE)
	
	ppmax <- function(m., max.m., lambda.){
		out <- dpois( m., lambda = lambda. ) 
		out[is.na(m.)] <- ppois( max.m., lambda = lambda., low = FALSE) 
		out
	}
	
	dt_pspec <- apsis_dat$PspecDT[[POP]]
	
	dt_pspec$pspec.est <- dt_pspec[[pr_pspec_var]]
	
	PR_PSPEC <- with(dt_pspec[,list(
		weight = ppmax( m.=m, max.m. = max.m, lambda. = 2 * N_total * maf_pop ), 
		pspec.est = pspec.est
	),], 
		sum(pspec.est*weight)/sum(weight)
	)
	
	if( N %% 100 == 0 ){
		tmp <- subset(apsis_dat$RsqList[[POP]][[ARR]][[round(N/100)+1]], Rsq>=rsq_thresh)[,list(
			weight = sum( ppmax( m.=m, max.m. = max.m, lambda. = 2 * N_total * maf_gwas ) * (
				p.Rsq.pspec*PR_PSPEC + p.Rsq.shared*(1-PR_PSPEC) )  
			) 
		), by=list(Rsq)]
	}else{
		n.a <- floor(N/100)
		n.b <- ceiling(N/100)
		w.a <- 1/abs(N/100 - n.a)/(1/abs(N/100 - n.a) + 1/abs(N/100 - n.b))

		tmp <- merge(
			x = subset(apsis_dat$RsqList[[POP]][[ARR]][[n.a + 1]], Rsq>=rsq_thresh)[,list(
				weight.a = w.a*sum(ppmax(m.=m,max.m.=max.m,lambda.=2*N_total*maf_gwas) * (	
					p.Rsq.pspec*PR_PSPEC + p.Rsq.shared*(1-PR_PSPEC) ) )
				), by=list(Rsq)], 
			y = subset(apsis_dat$RsqList[[POP]][[ARR]][[n.b + 1]], Rsq>=rsq_thresh)[,list(
				weight.b = (1-w.a)*sum(ppmax(m.=m,max.m.=max.m,lambda.=2*N_total*maf_gwas) * (
				p.Rsq.pspec*PR_PSPEC + p.Rsq.shared*(1-PR_PSPEC) ) )
		), by=list(Rsq)], 
			by = 'Rsq', all = TRUE)[,list(
			weight = sum(c(weight.a,weight.b),na.rm=TRUE)
		),by=Rsq]
	}
	
	power_out <- sapply(M, function(n.imp){
			with(tmp,
				sum( weight*pnorm( sqrt(n.imp*Rsq + N)*eta - c_a ))
			)
	})
	
	mean_rsq <- with(tmp,sum( weight*Rsq ))
	
	eff_type <- with(disease_model$model,paste(ifelse(is.na(rr),'LThres','RR')))
	if(disease_model$model$type != "") eff_type <- paste0(eff_type, ',', disease_model$model$type )
	
	data.table(
		# 'Model' = disease_model$model_text,
		'MAF' = disease_model$model$maf,
		'EffSize' = with(disease_model$model,ifelse(is.na(rr),beta,rr)),
		'Type' = eff_type,
		'Prev' = disease_model$model$prev,
		'Population' = POP_NAME,
		'Array' = ARR_NAME,
		'N_Sequenced' = N,
		'N_Imputed' = M,
		'Cost' = N*SEQ_COST + M*ARR_COST,
		'Power' = power_out,
		'E_Rsq' = mean_rsq
	)
}
