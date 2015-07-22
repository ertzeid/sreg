


`%ni%` <- Negate(`%in%`) 


logit = function(x){1/(1+exp(-x))}
probit = function(x){pnorm(x)}


#Define posterior simulation function
do2SGAM = function(fs.list,znames,m,N,n,kn=KK,fam='binom',d=data){
#fs.list = list(fs.BC,fs.share,fs.mean,fs.sd)
#m = mobs1
#z=c('res.BC14','res.share','res.mean','res.sd')
#znames=z
#N=25
#n=100
#kn=KK
#fam='binom'
#d<-data
	set.seed(413)
	B = m$coef
	library(MASS)
	for (i in 1:N){
		nd<-d
		for (j in 1:length(fs.list)){
			fs<-fs.list[[j]]
			brep<-mvrnorm(1,mu = coefficients(fs),Sigma =  vcov(fs))
			fs$coefficients<-brep
			zeta<-nd[,colnames(nd) == colnames(fs$model)[1]] -predict(fs,newdata = nd,type='response')
			nd[,colnames(nd) == znames[j]]<-zeta	
		}
		if (fam=="binom"){
			mssk =  gam(m$formula
				,knots=kn,family=binomial,data=nd,method='REML')
		}
		if (fam=="nb"){
			mssk =  gam(m$formula
				,knots=kn,family=nb(),data=nd,method='REML')
		}
		plot(mssk,pages=1,scheme=2)
		print(summary(mssk))
		Bs = mvrnorm(n,mu = mssk$coef, Sigma = mssk$Vp)
		B = cbind(B,t(Bs))
		print(paste(i,'of',N))
	}
	mf = m
	mf$Vp<-cov(t(B))
	return(mf)
}
do2SGAM = function(fs.list,znames,islog = NULL,m,N,n,kn=KK,fam='binom',d=data,plot.me=TRUE,link='logit'){
#fs.list = list(fs.BC,fs.share,fs.mean,fs.sd)
#m = mobs1
#z=c('res.BC14','res.share','res.mean','res.sd')
#znames=z
#N=25
#n=100
#kn=KK
#fam='binom'
#d<-data
	if (!is.null(islog)){
		if (length(fs.list) != length(islog)){stop('Fix the "islog" argument')}
	}
	set.seed(413)
	B = m$coef
	library(MASS)
	for (i in 1:N){
		nd<-d
		for (j in 1:length(fs.list)){
			fs<-fs.list[[j]]
			brep<-mvrnorm(1,mu = coefficients(fs),Sigma =  vcov(fs))
			fs$coefficients<-brep
			if (islog[j] == 0){
				zeta<-nd[,colnames(nd) == colnames(fs$model)[1]] -predict(fs,newdata = nd,type='response')
			}else{
				varname = substr(colnames(fs$model)[1],start = 5,stop = nchar(colnames(fs$model)[1])-1)
				zeta<-nd[,colnames(nd) == varname] -exp(predict(fs,newdata = nd,type='response'))
			}
			nd[,colnames(nd) == znames[j]]<-zeta	
		}
		if (fam=="binom"){
			mssk =  gam(m$formula
				,knots=kn,family=binomial(link=link),data=nd,method='REML')
		}
		if (fam=="gaussian"){
			mssk =  gam(m$formula
				,knots=kn,data=nd,method='REML')
		}
		if (fam=="nb"){
			mssk =  gam(m$formula
				,knots=kn,family=nb(),data=nd,method='REML')
		}
		if (plot.me == TRUE){
			plot(mssk,pages=1,scheme=2)
		}
		print(summary(mssk))
		Bs = mvrnorm(n,mu = mssk$coef, Sigma = mssk$Vp)
		B = cbind(B,t(Bs))
		print(paste(i,'of',N))
	}
	mf = m
	mf$Vp<-cov(t(B))
	return(mf)
}

library(plyr)
rep.row <- function(r, n){
  colwise(function(x) rep(x, n))(r)
}

NOIVASF = function(model, term, cond){
#model = ss0
#term = 'TP.link.imp.mean'
#cond = list(TP=0,FT=0,price_lr14=1250)

	m<-model
	d<-m$model
	x<-d[,colnames(d) == term]
	if (length(cond)>0){
		for (i in 1:length(cond)){
			d[,colnames(d) == names(cond)[i]]<-cond[[i]]
		}
	}
	means = as.data.frame(t(colMeans(d)))
	xgrid = with(d, data.frame(x = seq(from=min(x),to=max(x),length.out=50)))
	xgrid$E<-xgrid$ciup<-xgrid$cidown<-NA
	for (i in 1:nrow(xgrid)){
		means[,colnames(means) == term]<-xgrid$x[i]
		p<-predict(m,newdata=means,se.fit=T)
		xgrid$E[i] = pnorm(p$fit)
		xgrid$ciup[i] = pnorm(p$fit+1.96*p$se.fit)
		xgrid$cidown[i] = pnorm(p$fit-1.96*p$se.fit)
	}
	return(xgrid)
}

vcovHC_GAM = function(m){
	bread = m$Vp/m$sig2
	X = predict(m,type='lpmatrix')
	e = residuals(m)
	omega <- function(residuals, diaghat, df) residuals^2 
#	meat = t(X)%*%e%*%t(e)%*%X
	meat = t(X)%*%diag(omega(e))%*%X
	V = bread%*%meat%*%bread
	N = nrow(m$model)
	K = sum(m$edf)
	return((N-1)/(N-K)*V)
}

makemeanframe = function(m,gridvar,gridvarval,resvars,cond=list()){#function to return a data frame of all covariates held to their means except the first stage residuals and the grid variable
	d = m$model
	rvs = matrix(rep(NA,length(resvars)*nrow(d)),nrow(d))
	for (i in 1:length(resvars)){
		rvs[,i]<-d[,colnames(d) == resvars[i]]
	}
	colnames(rvs)<-resvars
	mm = colMeans(d)
	mm = rep.row(as.data.frame(t(mm)),nrow(d))
	mm[,colnames(mm) == gridvar]<-gridvarval
	mm = mm[,colnames(mm) %ni% resvars]
	mm = cbind(mm,rvs)
	if (length(cond)>0){
		for (i in 1:length(cond)){
			mm[,colnames(mm) == names(cond)[i]]<-cond[[i]]
		}
	}
	lpm = predict(m, type='lpmatrix',newdata=mm)
	return(lpm)
}

getASF = function(m,var,resvars,gridlength = 50,cond=list(),B = coef(m),X = NULL){#function to calculate ASF from a 2-stage (CF) GAM
	v = m$model[,colnames(m$model) == var]
	asf <- with(m$model, cbind(seq(from=min(v,na.rm=TRUE),to=max(v,na.rm=TRUE),length.out = gridlength),NA))
	for(i in 1:dim(asf)[1]){
		if (is.null(X)){
			dat <- makemeanframe(m,var,asf[i,1],resvars,cond=cond)
		} else{
			dat<-X[[i]]
		}
		asf[i,2]<-mean(pnorm(dat%*%B))
	}
	return(asf)
}

simbootASF = function(m,n=1000,var,resvars,gridlength = 50,cond=list(),verbose=F){
#n=1000
#m=ss.sim
#var = 'FT'
#resvars = c('res.tpmean','res.tpsd')
#gridlength=2
#cond=list(TP=0,price_lr14=1250)
	library(MASS)
	Breps = mvrnorm(n,mu = m$coef, Sigma = m$Vp)
	v = m$model[,colnames(m$model) == var]
	xgrid = seq(from=min(v,na.rm=TRUE),to=max(v,na.rm=TRUE),length.out = gridlength)
	ASFmat = matrix(rep(NA,n*gridlength),gridlength)
	X = list()
	for (i in 1:gridlength){
		X[[i]] = makemeanframe(m,var,xgrid[i],resvars,cond=cond)
	}
	if (verbose == T){print('made meanframes')}
	for (i in 1:n){
		ASFmat[,i] = getASF(m,var,resvars,gridlength,cond=cond,B = Breps[i,],X)[,2]
		if (verbose == T){print(i)}
	}
	E = rowMeans(ASFmat)
	ciup = apply(ASFmat,1,quantile,probs = .975)
	cidown = apply(ASFmat,1,quantile,probs = .025)
	return(list(x = xgrid,E = E,ciup = ciup,cidown = cidown))
}

vcovHC_GAM = function(m){
	bread = m$Vp/m$sig2
	X = predict(m,type='lpmatrix')
	e = residuals(m)
	omega <- function(residuals, diaghat, df) residuals^2 
#	meat = t(X)%*%e%*%t(e)%*%X
	meat = t(X)%*%diag(omega(e))%*%X
	V = bread%*%meat%*%bread
	N = nrow(m$model)
	K = sum(m$edf)
	return((N-1)/(N-K)*V)
}




#####################
#Special regressors
sreg = function(y,ex,en,I, V, data,bw = "nrd0", qupchop = 1, qdownchop = 0, absupchop = Inf, absdownchop = -Inf){
	D<-data
	y= D[,colnames(D)==y]
	S = paste(ex,en,I)
	X = paste(ex,en)
	Z = paste(ex,I)

	m = lm(as.formula(paste("V~",S)),data=D)
	D$Uhat = V - predict(m, newdata = D)

	f = density(D$Uhat,na.rm=TRUE,bw=bw)
	getfhat <- with(f, approxfun(x, y, rule=1))
	D$fhat = getfhat(D$Uhat)

	D$That = with(D, (y - (V>=0))/fhat)

	D = subset(D, That >=quantile(D$That,probs = qdownchop,na.rm=T) & That <=quantile(D$That,probs = qupchop,na.rm=T))
	D = subset(D, That >=absdownchop & That <=absupchop)

	library(AER)
	m = ivreg(as.formula(paste('That~',X,'|',Z)),data=D)
	
	return(list(m = m, y=y,V = V, Uhat = D$Uhat, f = f, fhat = D$fhat ,That = D$That, D = D))
}

sregGAM = function(y, ex, en, I, resvars, V,
                   fs.list = list(), islog=NULL, znames, data,
                   bw = "nrd0", qupchop = 1, qdownchop = 0,
                   absupchop = Inf, absdownchop = -Inf,
                   sim=TRUE, verbose = FALSE,allcrossprod = F,Vname = NULL){
#en = "+s(TP.link.imp.mean,k=K)+s(TP.link.imp.sd,k=K)"
#I = "+s(TPmz_y.mean,k=K)+s(TPmz_y.sd)"
#resvars = "+s(res.tpmean,k=K)+s(res.tpsd,k=K)"
#V='V'
#y='BC14'
#znames = c('res.tpmean','res.tpsd')
#islog=c(0,1)
#V = V
# fs.list = list(sr.fs1.tpmean,sr.fs1.tpsd)
#znames = znames
#data=srdat
#sim=T
	library(mgcv)
	library(MASS)
	D<-data
	D$y = D[,colnames(D)==y]
	S = paste(ex, en, I)
	X = paste(ex, en)
	Z = paste(ex, I)
	if (allcrossprod == T){#Augments S with all cross-products
		m = gam(as.formula(paste(V, "~", S)), data=D, knots=KK, method='REML',fit=F)
		varnames = names(m$var.summary)
		varnames=varnames[varnames != 'Vname']
		counter = 0
		for (i in 1:(length(varnames)-1)){
			v1 = D[,colnames(D)== varnames[i]]
			for (j in (i+1):length(varnames)){
				v2 = D[,colnames(D) == varnames[j]]
				counter = counter+1
				D$vnew = v1*v2
				colnames(D)[colnames(D)=='vnew']<-paste('CP_',counter,sep='')
				S = paste(S,paste('+CP_',counter,sep=''))
			}
		}
	}
	m = gam(as.formula(paste(V, "~", S)), data=D, knots=KK, method='REML')
	if (verbose == TRUE) { print('fit V = f(S)') }
	m$Vp<-vcovHC_GAM(m)#adjust for heteroskedasticity
	D$Uhat = D$V - predict(m, newdata = D)
	f = density(D$Uhat,na.rm=TRUE,bw=bw)
	getfhat <- with(f, approxfun(x, y, rule=1))
	D$fhat = getfhat(D$Uhat)
	D$That = with(D, (D$y - (V>=0))/fhat)

	D = subset(D, That >=quantile(D$That,probs = qdownchop,na.rm=T) & That <=quantile(D$That,probs = qupchop,na.rm=T))
	D = subset(D, That >=absdownchop & That <=absupchop)

	ss.sreg = gam(as.formula(paste('That~',ex,en,resvars)), family=gaussian,method='REML',data=D,knots=KK)
	ss.sreg$Vp<-vcovHC_GAM(ss.sreg)#adjust for heteroskedasticity
	m<-ss.sreg
	if (sim==TRUE){
		ss.sreg.sim = do2SGAM(fs.list, islog = islog, znames = znames, m=ss.sreg, N=25, n=1000, fam='gaussian', d=D)
		m<-ss.sreg.sim
	}
	return(list(m = m, y=D$y, V = D$V, Uhat = D$Uhat, f = f, fhat = D$fhat ,That = D$That, D = D,znames = znames))
}

getAIF = function(sreg.obj,var,resvars,gridlength = 50,cond=list(),B = NULL,X = NULL,getforV = F,method='bGAM',loess.edf,functogetV = NULL,invfuncV = NULL,Vvar = 'price_lr14'){#function to calculate AIF from a 2-stage (CF) SREG
#cond = list(TP=0,FT=0)
#var = 'V'
	m = sreg.obj$m
	if (is.null(B)){B = coef(m)}
	lpm = predict(m,type='lpmatrix')
	XB = lpm%*%B+sreg.obj$V
	if (method == 'bGAM'){l = gam(sreg.obj$y~s(XB,k=20),knots = list(XB = seq(from=min(predict(m)) ,to = max(predict(m)) ,length.out=20)),family=binomial,method='REML')}
	if (method == 'loess'){l = loess(sreg.obj$y~XB,enp.target=loess.edf)}
	if (Vvar %in% names(cond)){
		V<-functogetV(cond[names(cond) == Vvar][[1]])
	}else{
		V = sreg.obj$V
	}
	v = sreg.obj$D[,colnames(sreg.obj$D) == var]
	aif <- with(m$model, cbind(seq(from=min(v,na.rm=TRUE),to=max(v,na.rm=TRUE),length.out = gridlength),NA))
	for(i in 1:dim(aif)[1]){
		if (is.null(X)){
			dat <- makemeanframe(m,var,aif[i,1],resvars,cond=cond)
		} else{
			dat<-X[[i]]
		}
		if (getforV == T){#Either add V (which gets averaged over) or add a value along its gradient
			XBi = dat%*%B+functogetV(aif[i,1])
			}else{
			XBi = dat%*%B+V
		}
		aif[i,2] <- predict(l, newdata=data.frame(XB=mean(XBi)),type='response')
	}
	return(list(aif=aif,index=l))
}

simbootAIF = function(sreg.obj,n=1000,var,resvars,gridlength = 50,cond=list(),verbose=F,Vvar = "price_lr14",functogetV = NULL, invfuncV = NULL,method = 'bGAM',loess.edf=5){
#sreg.obj = sreg1.sim
#n=10
#var = 'TP.link.imp.mean'
#resvars=resvars
#gridlength = 50
#cond=list(TP=0,FT=0,price_lr14 = log(100))
#verbose=T
#Vvar = 'price_lr14'
#functogetV = function(x){-(x-mean(srdat$price_lr14,na.rm=TRUE))}
#invfuncV = function(V){mean(srdat$price_lr14,na.rm=TRUE)-V}
#method = 'bGAM'
#loess.edf=6
#rm(sreg.obj,n,var,Vvar,resvars,gridlength,cond,verbose,functogetV,invfuncV,getforV)
	m<-sreg.obj$m
	library(MASS)
	Breps = mvrnorm(n,mu = m$coef, Sigma = m$Vp)
	V = sreg.obj$V
	if (Vvar %in% names(cond)){
		if (is.null(functogetV)){stop("Dogg.  You saying that V ain't even transformed?")}
		V<-functogetV(cond[names(cond) == Vvar][[1]])
	}
	if (var == Vvar){
		if (is.null(functogetV)){stop("Dogg.  You saying that V ain't even transformed?")}
		if (is.null(invfuncV)){stop("Dogg.  You ain't tell me how to go from V to whatever.")}
		getforV<-TRUE
	}else{
		getforV<-FALSE
	}
	v = sreg.obj$D[,colnames(sreg.obj$D) == var]

	xgrid = seq(from=min(v,na.rm=TRUE),to=max(v,na.rm=TRUE),length.out = gridlength)
	AIFmat = matrix(rep(NA,n*gridlength),gridlength)
	X = list()
	if (var == Vvar){#The mean frames are identical when making a grid over the special regressor, because the value of the special regressor along the grid is simply added into the index function in the "getAIF" function -- this if statement saves time.
		x = makemeanframe(m,var,xgrid[1],resvars = sreg.obj$znames,cond=cond)		
		for (i in 1:gridlength){X[[i]]<-x}
	}else{	
		for (i in 1:gridlength){
			X[[i]] = makemeanframe(m,var,xgrid[i],resvars = sreg.obj$znames,cond=cond)
		}
	}
	if (verbose == T){print('made meanframes')}
	p<-beta<-list() #for returning index functions
	for (i in 1:n){#Get AIFs from n (parametric) bootstrap samples

		AIF = getAIF(sreg.obj,var,resvars,gridlength,cond=cond,B = Breps[i,]
			,X=X,getforV = getforV,functogetV = functogetV,invfuncV = invfuncV,Vvar=Vvar,method = method,loess.edf = loess.edf)
		AIFmat[,i] = AIF$aif[,2]
		if (method == 'bGAM'){beta[[i]]<-AIF$index$coef}
		if (method == 'loess'){p[[i]]<-predict(AIF$index, newdata = data.frame(XB = qnorm(seq(from=.0001,to=.9999,length.out=1001))))}
		if (verbose == T){print(i)}
	}
	if (method == 'bGAM'){
		bmat = matrix(unlist(beta),length(beta[[1]]))
		bhat<-rowMeans(bmat)
		Vp<-cov(t(bmat))
		index<-AIF$index
		index$coefficients<-bhat
		index$Vp<-Vp
		xgrid_index = data.frame(XB = seq(from=min(predict(m)),to=max(predict(m)),length.out=1001))
		p = predict(index,newdata = xgrid_index,se.fit=TRUE)
		xgrid_index$y = pnorm(p$fit)
		xgrid_index$ciup = pnorm(p$fit+p$se.fit*1.96)
		xgrid_index$cidown = pnorm(p$fit-p$se.fit*1.96)
		indexDat = list(xgrid = xgrid_index,index = index)
	}
	if (method == 'loess'){
		xgrid_index = data.frame(XB = seq(from=min(predict(m)),to=max(predict(m)),length.out=1001))
		pmat = matrix(unlist(p),length(p[[1]]))
		xgrid_index$y = rowMeans(pmat,na.rm=TRUE)
		xgrid_index$ciup = apply(pmat,1,quantile,probs=.975,na.rm=TRUE)
		xgrid_index$cidown = apply(pmat,1,quantile,probs=.025,na.rm=TRUE)
		indexDat = list(xgrid = xgrid_index)
	}
	E = rowMeans(AIFmat,na.rm=TRUE)
	ciup = apply(AIFmat,1,quantile,probs = .975,na.rm=TRUE)
	cidown = apply(AIFmat,1,quantile,probs = .025,na.rm=TRUE)
	return(list(x = xgrid,E = E,ciup = ciup,cidown = cidown, indexDat = indexDat))
}


pp = function(model, term){
	xpos = par("usr")[1]+(par("usr")[2]-par("usr")[1])*.3
	ypos = par("usr")[4]-(par("usr")[4]-par("usr")[3])*.15
	pval = round(summary(model)$s.table[term,4],2)
	if (pval<.01) { pval<-'.01';sym<-'<'}else{sym='='}
	string = paste('p',sym,pval,sep='')
	text(xpos,ypos,string,col='blue')
}


makearrows = function(m, v, y, x1 = mean(m$model[,colnames(m$model)==v]), x2 = x1+sd(m$model[,colnames(m$model)==v]),yoff=1,xoff=0){
	segments(x1, y,x2, y, col='red')
	x = m$model
	x[colnames(x) == v]<-x1
	l = predict(m,newdata=x)
	x[colnames(x) == v]<-x2
	h = predict(m,newdata=x)
	or = round(mean(exp(h-l)),2)
	text(mean(c(x1,x2))+xoff,y+yoff, or)
}

