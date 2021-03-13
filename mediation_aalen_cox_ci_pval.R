## method=="Aalen" -> output Hazard Difference
## method=="Cox"   -> output Hazard Ratio

mediation_ci1 <- function(lambda.s, lambda.g, covar11, covar12,
	covar22, alpha.s, var_alpha, G=10^4, method){
	require(mvtnorm)
	Omega <- matrix(c(covar11,covar12,covar12,covar22),nrow=2)
	IE <- rep(0,G); DE <- rep(0,G); TE <- rep(0,G); Q <- rep(0,G)

	set.seed(137)
	lambda <- rmvnorm(G, mean = c(lambda.s, lambda.g), sigma = Omega)
	alpha <- rnorm(G, mean=alpha.s, sd=sqrt(var_alpha))
	DE <- lambda[,1]
	IE <- lambda[,2] * alpha
	TE <- IE + DE

	DE.obs <- lambda.s
	IE.obs <- lambda.g * alpha.s
	TE.obs <- DE.obs+IE.obs
	pval.DE<-2*min(mean((DE-mean(DE))>DE.obs), mean((DE-mean(DE))<DE.obs))
	pval.IE<-2*min(mean((IE-mean(IE))>IE.obs), mean((IE-mean(IE))<IE.obs))
	pval.TE<-2*min(mean((TE-mean(TE))>TE.obs), mean((TE-mean(TE))<TE.obs))

	if (method=="Cox") {DE=exp(DE); IE=exp(IE); TE=exp(TE)}
	print("DE:")
	print(ifelse(method=="Aalen", DE.obs, exp(DE.obs)))
	print(quantile(DE, c(0.025, 0.975)))
	print(paste("pval_DE=", pval.DE))
	print("IE:")
	print(ifelse(method=="Aalen", IE.obs, exp(IE.obs)))
	print(quantile(IE, c(0.025, 0.975)))
	print(paste("pval_IE=", pval.IE))
	print("TE:")
	print(ifelse(method=="Aalen", TE.obs, exp(TE.obs)))
	print(quantile(TE, c(0.025, 0.975)))
	print(paste("pval_TE=", pval.TE))
}

mediation_ci2 <- function(lambda.s, lambda.m, lambda.g, Sigma.lambda,
	alpha.s, alpha.m, Sigma.alpha, delta.s, Sigma.delta,
	G=10^4, method){
	require(mvtnorm)
	SY <- rep(0,G); SGY <- rep(0,G); SMY <- rep(0,G); TE <- rep(0,G)
	set.seed(137)
	lambda <- rmvnorm(G, mean = c(lambda.s, lambda.m, lambda.g),
		sigma = Sigma.lambda)
	alpha <- rmvnorm(G, mean = c(alpha.s, alpha.m),
		sigma = Sigma.alpha)
	delta <- rnorm(G, mean=delta.s, sd=sqrt(Sigma.delta))
	SY <- lambda[,1]
	SGY <- lambda[,3] * alpha[,1]
	SMY <- (lambda[,2] + lambda[,3]*alpha[,2])*delta
	TE <- SY+SGY+SMY

	SY.obs <- lambda.s
	SGY.obs <- lambda.g * alpha.s
	SMY.obs <- (lambda.m + lambda.g*alpha.m)*delta.s
	TE.obs <- SY.obs+SGY.obs+SMY.obs
	pval.SY<-2*min(mean((SY-mean(SY))>SY.obs), mean((SY-mean(SY))<SY.obs))
	pval.SGY<-2*min(mean((SGY-mean(SGY))>SGY.obs), mean((SGY-mean(SGY))<SGY.obs))
	pval.SMY<-2*min(mean((SMY-mean(SMY))>SMY.obs), mean((SMY-mean(SMY))<SMY.obs))
	pval.TE<-2*min(mean((TE-mean(TE))>TE.obs), mean((TE-mean(TE))<TE.obs))

	if (method=="Cox") {SY=exp(SY); SGY=exp(SGY); SMY=exp(SMY); TE=exp(TE)}
	print("SY:")
	print(ifelse(method=="Aalen", SY.obs, exp(SY.obs)))
	print(quantile(SY, c(0.025, 0.975)))
	print(paste("pval_SY=", pval.SY))
	print("SGY:")
	print(ifelse(method=="Aalen", SGY.obs, exp(SGY.obs)))
	print(quantile(SGY, c(0.025, 0.975)))
	print(paste("pval_SGY=", pval.SGY))
	print("SMY:")
	print(ifelse(method=="Aalen", SMY.obs, exp(SMY.obs)))
	print(quantile(SMY, c(0.025, 0.975)))
	print(paste("pval_SMY=", pval.SMY))
	print("TE:")
	print(ifelse(method=="Aalen", TE.obs, exp(TE.obs)))
	print(quantile(TE, c(0.025, 0.975)))
	print(paste("pval_TE=", pval.TE))
}
