KS.prelim <- function(Y, X=NULL, id=NULL, out_type = "D", B=1000, model = "gaussian"){
  ##Preliminary
  Y = as.matrix(Y);n = nrow(Y)

  if(length(X)!=0){X0 = svd(as.matrix(X))$u}else{X0 = NULL}
  X0 = cbind(rep(1,n),X0)

  # if(out_type=="C"){nullglm = glm(Y~0+X0,family=gaussian)}
  # if(out_type=="D"){nullglm = glm(Y~0+X0,family=binomial)}

  nullglm = glm(Y~0+X0, family = model)

  if (length(id)==0){id = 1:n}

  mu = nullglm$fitted.values;Y.res = Y-mu
  #permute the residuals
  index = sapply(1:B,function(x)sample(1:length(Y)));temp.Y.res = Y.res[as.vector(index)]
  re.Y.res = matrix(temp.Y.res,length(Y),B)

  #prepare invserse matrix for covariates
  if(out_type=='D'){v = mu*(1-mu)}else{v = rep(as.numeric(var(Y.res)),length(Y))}
  inv.X0 = solve(t(X0)%*%(v*X0))

  #prepare the preliminary features
  result.prelim = list(Y=Y,id=id,n=n,X0=X0,nullglm=nullglm,out_type=out_type,re.Y.res=re.Y.res,inv.X0=inv.X0)
  return(result.prelim)
}

Get.p <- function(X,result.prelim){
  X = as.matrix(X)
  mu = result.prelim$nullglm$fitted.values;Y.res = result.prelim$Y-mu
  outcome = result.prelim$out_type
  if(outcome=="D"){
    p = ScoreTest_SPA(t(X),result.prelim$Y,result.prelim$X0,method=c("fastSPA"),minmac=-Inf)$p.value
  }else{
    v = rep(as.numeric(var(Y.res)),nrow(X))
    p = pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
  }
  return(as.matrix(p))
}

Get.p.base <- function(X,result.prelim){
  X = as.matrix(X)
  mu = result.prelim$nullglm$fitted.values;Y.res = result.prelim$Y-mu
  outcome = result.prelim$out_type
  if(outcome=="D"){v = mu*(1-mu)}else{v = rep(as.numeric(var(Y.res)),nrow(X))}
  p = pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
  p[is.na(p)] = NA
  return(p)
}

Get.p.moment <- function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
  re.mean = apply(re.Q,2,mean)
  re.variance = apply(re.Q,2,var)
  re.kurtosis = apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
  re.df = (re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
  re.p = t(pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F))
  #re.p[re.p==1] = 0.99
  return(re.p)
}

SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  r  =  length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  out  =  .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
  out$res  =  1 - out$res
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
}

Get_Liu_PVal.MOD.Lambda <- function(Q.all, lambda, log.p=FALSE){
  param = Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm = (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 = Q.Norm * param$sigmaX + param$muX
  p.value =  pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  return(p.value)
}

Get_Liu_Params_Mod_Lambda <- function(lambda){
  ## Helper function for getting the parameters for the null approximation

  c1 = rep(0,4)
  for(i in 1:4){
    c1[i] = sum(lambda^i)
  }

  muQ = c1[1]
  sigmaQ = sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1 = sqrt(8)*s1
  beta2 = 12*s2
  type1 = 0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1 = 1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX  = l+d
  sigmaX = sqrt(2) *a

  re = list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get.p.SKAT <- function(score,re.score,K,window.matrix,weight,result.prelim){

  mu = result.prelim$nullglm$fitted.values;Y.res = result.prelim$Y-mu
  X0 = result.prelim$X0;outcome = result.prelim$out_type

  Q = as.vector(t(score^2)%*%(weight*window.matrix)^2)
  K.temp = weight*t(weight*K)

  #fast implementation by resampling based moment matching
  p0 = Get.p.moment(as.vector(t(score^2)%*%(weight*window.matrix)^2),re.score^2%*%(weight*window.matrix)^2)
  p = p0#c()
  for(i in which(p0<0.01 |p0>=1)){
    #print(i)
    temp = K.temp[window.matrix[,i]!=0,window.matrix[,i]!=0]
    if(sum(temp^2)==0){p[i] = NA;next}

    #proc.time()
    lambda=eigen(temp,symmetric=T,only.values=T)$values
    temp.p = SKAT_davies(Q[i],lambda,acc=10^(-6))$Qq
    #proc.time()

    #proc.time()
    #temp.p = Get_Liu_PVal.MOD(Q[i],temp*2)$p.value
    #proc.time()

    #temp.p = davies(Q[i],lambda,acc=10^(-6))$Qq
    if(temp.p > 1 || temp.p <= 0 ){
      temp.p = Get_Liu_PVal.MOD.Lambda(Q[i],lambda)
    }
    p[i] = temp.p
  }
  return(as.matrix(p))
}

Get.cauchy <- function(p){
  p[p>0.99] = 0.99
  is.small = (p<1e-16) & !is.na(p)
  is.regular = (p>=1e-16) & !is.na(p)
  temp = rep(NA,length(p))
  temp[is.small] = 1/p[is.small]/pi
  temp[is.regular] = as.numeric(tan((0.5-p[is.regular])*pi))

  cct.stat = mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}

Get.cauchy.scan <- function(p,window.matrix){
  p[p>0.99] = 0.99
  is.small = (p<1e-16)
  temp = rep(0,length(p))
  temp[is.small] = 1/p[is.small]/pi
  temp[!is.small] = as.numeric(tan((0.5-p[!is.small])*pi))
  #window.matrix.MAC10 = (MAC>=10)*window.matrix0

  cct.stat = as.numeric(t(temp)%*%window.matrix/apply(window.matrix,2,sum))
  #cct.stat = as.numeric(t(temp)%*%window.matrix.MAC10/apply(window.matrix.MAC10,2,sum))
  is.large = cct.stat>1e+15 & !is.na(cct.stat)
  is.regular = cct.stat<=1e+15 & !is.na(cct.stat)
  pval = rep(NA,length(cct.stat))
  pval[is.large] = (1/cct.stat[is.large])/pi
  pval[is.regular] = 1-pcauchy(cct.stat[is.regular])
  return(pval)
}

KS_test <- function(G, MAF, MAC, result.prelim, window.matrix, weight.matrix){
  # Burden test
  p.burden = matrix(NA,ncol(window.matrix),ncol(weight.matrix))
  for (k in 1:ncol(weight.matrix)){
    temp.window.matrix = weight.matrix[,k]*window.matrix
    X = as.matrix(G%*%temp.window.matrix)
    p.burden[,k] = Get.p.base(X,result.prelim)
  }

  # Dispersion test
  mu = result.prelim$nullglm$fitted.values
  Y.res = result.prelim$Y-mu
  re.Y.res = result.prelim$re.Y.res
  X0 = result.prelim$X0
  outcome = result.prelim$out_type
  if(outcome=="D"){v = mu*(1-mu)}else{v = rep(as.numeric(var(Y.res)),nrow(G))}
  A = t(G)%*%(v*G)
  B = t(G)%*%(v*X0)
  C = solve(t(X0)%*%(v*X0))
  K = A-B%*%C%*%t(B)
  score = t(G)%*%Y.res
  re.score = t(t(G)%*%re.Y.res)
  p.dispersion = matrix(NA,ncol(window.matrix),ncol(weight.matrix))
  for (k in 1:ncol(weight.matrix)){
    p.dispersion[,k] = Get.p.SKAT(score,re.score,K,window.matrix,weight=(MAC>=5)*weight.matrix[,k],result.prelim)
  }

  # single variant test
  p.single = Get.p(G,result.prelim)

  # combination of single rare/common variants
  p.V1 = Get.cauchy.scan(p.single,(MAC>=5 & MAF<0.01)*window.matrix)
  p.V2 = Get.cauchy.scan(p.single,(MAF>=0.01)*window.matrix)

  p.individual = cbind(p.burden,p.dispersion,p.V1,p.V2);
  colnames(p.individual) = c(paste0("burden_",colnames(weight.matrix)),paste0("dispersion_",colnames(weight.matrix)),"singleCauchy_MAF<0.01&MAC>=5","singleCauchy_MAF>=0.01")

  # Cauchy combination of all tests
  p.KS = as.matrix(apply(p.individual,1,Get.cauchy))
  # Cauchy combination of all tests derived from common variants
  test.common = grep("MAF>=0.01",colnames(p.individual))
  p.KS.common = as.matrix(apply(p.individual[,test.common,drop=FALSE],1,Get.cauchy))
  # Cauchy combination of all tests derived from rare variants
  p.KS.rare = as.matrix(apply(p.individual[,-test.common,drop=FALSE],1,Get.cauchy))

  return(list(p.KS=p.KS, p.KS.common=p.KS.common,p.KS.rare=p.KS.rare,p.individual=p.individual))
}


