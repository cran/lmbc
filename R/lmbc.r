lmbc <- function(data,up,down,power)
{
	data <- expData2nt(g1_part, up, down) 
	numGenes <- length(unique(data$index))      # number of genes in dataset

	createPosition <- function(numGenes,up,down){
		position = numeric(0)

		for(i in up:1)
		{
			position = append(position,rep(i,3))
		}
		position = append(position,rep(0,3))
		for(i in 1:(down-1))
		{
			position = append(position,rep(i,3))
		}
		position + 1
	}

	lasso.adapt.bic2 <- function(x,y,z,position){

		# adaptive lasso from lars with BIC stopping rule 
		# this one uses the "known variance" version of BIC with RSS/(full model mse)
		# must use a recent version of R so that normalize=FALSE can be used in lars

		require(lars)
		ok<-complete.cases(x,y)
		x<-x[ok,]                            # get rid of na's
		y<-y[ok]                             # since regsubsets can't handle na's
		m<-ncol(x)
		n<-nrow(x)
		x<-as.matrix(x)                      # in case x is not a matrix

		#  standardize variables like lars does 
		one <- rep(1, n)
		meanx <- drop(one %*% x)/n
		xc <- scale(x, meanx, FALSE)         # first subtracts mean
		normx <- sqrt(drop(one %*% (xc^2)))
		names(normx) <- NULL
		xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)
		
		out.ls=lm(y~xs)                      # ols fit on standardized
		beta.ols=out.ls$coeff[2:(m+1)]       # ols except for intercept
		w=(position)^z
		xs=scale(xs,center=FALSE,scale=w)  # xs times the weights
		object=lars(xs,y,type="lasso",normalize=FALSE)

		# get min BIC
		# bic=log(n)*object$df+n*log(as.vector(object$RSS)/n)   # rss/n version
		sig2f=summary(out.ls)$sigma^2        # full model mse
		bic2=log(n)*object$df+as.vector(object$RSS)/sig2f       # Cp version
		step.bic2=which.min(bic2)            # step with min BIC

		fit=predict.lars(object,xs,s=step.bic2,type="fit",mode="step")$fit
		coeff=predict.lars(object,xs,s=step.bic2,type="coef",mode="step")$coefficients
		coeff=coeff*w/normx                  # get back in right scale
		st=sum(coeff !=0)                    # number nonzero
		mse=sum((y-fit)^2)/(n-st-1)          # 1 for the intercept
		
		# this next line just finds the variable id of coeff. not equal 0
		if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind<-0
		return(list(fit=fit,st=st,mse=mse,x.ind=x.ind,coeff=coeff,object=object,
					bic2=bic2,step.bic2=step.bic2))
					
	}

	y <- log(data$count+1) 
	index <- data$index
	numGenes <- length(unique(index))
	data_back <- data[,-c(1,2)] 
	options(contrasts=c("contr.sum","contr.poly"))
	m <- lm(y~factor(index))
	data_mean <- model.matrix(m)
	data_total <- cbind(data_mean,data_back)
	data_total <- data_total[-1]

	position <- createPosition(numGenes, up, down)

	geneMean <- rep(0,100)
	for(i in 1:100)
	{
		geneMean[i] <- mean(y[data$index==i])
	}

	y_nomean <- (y - geneMean[data$index])


	data_back_single <- data_total[,-c(1:99,343:1062)]

	fit.lasso <- lasso.adapt.bic2(data_back_single,y_nomean,power,position)

	aa <- fit.lasso$coef != 0
	bb <- c(fit.lasso$coef[-c(1)] != 0, FALSE)
	cc <- c(fit.lasso$coef[-c(12)] != 0, FALSE, FALSE)
	dd <- rep(0, length(aa))

	for(i in 1:length(aa))
	{
		dd[i] <- aa[i]+bb[i]+cc[i]
	}

	dd <- dd == 3

	f <- min((which(dd==1)+2)[which((which(dd==1)+2) %% 3 ==0 )]-2)
	g <- max((which(dd==1)+2)[which((which(dd==1)+2) %% 3 ==0 )])

	ff <- 3*(up+down)+3*(f-1)+1
	gg <- 3*(up+down)+9*((g/3)-1)

	lengthUP <- 40-(f-1)/3
	lengthDOWN <- (g/3)-41
	seqLen <- lengthUP + lengthDOWN
	
	f <- f+99 
	g <- g+99
	ff <- ff+99 
	gg <- gg+99
	
	fit <- lm(y~.,data=data_total[c(1:99,f:g,ff:gg)])
	fit0 <- lm(y~.,data=as.data.frame(data_mean[,-1]))
	p1 <- predict.lm(fit,data_total[c(1:99,f:g,ff:gg)])
	p0 <- predict.lm(fit0,as.data.frame(data_mean[,-1]))
	r2 <- 1-(sum((y-p1)^2))/(sum((y-p0)^2))
	logLik <- logLik(fit)

	list(seqLen=seqLen,lengthUP=lengthUP,lengthDOWN=lengthDOWN,r2=r2,logLik=logLik)

}
