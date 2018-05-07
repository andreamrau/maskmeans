
#-------------------------------------------------------------------
## OLD COMMANDS: KEEPING THESE FOR BACKWARDS COMPATIBILITY FOR NOW
#-------------------------------------------------------------------


############################################################################################################
############################################################################################################
##   Splitting de clusters inspirée du multiview  -   poids par cluster 
## Entrées : 
#     probapost.init : une classif floue
#     X : les données (le 1er bloc correspond aux données principales utilisées pour obtenir classif)
#     mv : la taille de chaque multiview
#     Kmax : le nombre max de classes
#     gamma : la puissance sur les poids w
#     delta : la puissance sur les \pi_{i,k}

splittingProbapostbis=function(X,mv,gamma=2,delta=2,Kmax,probapost.init){
  #mettre ici un stop si delta <=1 ....
  
  ref=c(0,cumsum(mv))
  CRIT<-NULL    #evolution critere general
  probapost=probapost.init
  Kinit=ncol(probapost.init)
  
  # initialisation des centres
  centers<-matrix(0,nrow=Kinit,ncol=ncol(X))
  for (k in 1:Kinit){
    centers[k,]<- apply(matrix(rep(probapost[,k]^delta,sum(mv)),nrow=nrow(X)) * X,2,sum) / sum(probapost[,k]^delta)
  }
  rownames(centers)<-1:Kinit
  
  # calcul des poids initiaux
  weights.sauv<-NULL
  w=weightprobapostbis(X, mv,centers,probapost,gamma,delta)
  weights.sauv[[1]] = w$weights
  
  #calcul des variances pondérées par classe   -> optimiser ces boucles for!!
  varclass=rep(0,Kinit)
  for (k in 1:Kinit){
    for (i in 1:nrow(X)){
    for (v in 1:length(mv)){
      dv=(ref[v]+1):ref[v+1]
      varclass[k] = varclass[k] + ( (w$weights[k,v]^gamma) * (probapost[i,k]^delta) * sum( (as.matrix(X)[i,dv] - centers[k,dv])^2))
      }
    }
  }
  # calcul du critere general 
  CRIT<-c(CRIT, sum(varclass))
  
  ###   début de la boucle de splitting
  ksplit.sauv<-NULL
  iter=1
  nbcluster=Kinit
  
  while(nbcluster < Kmax ){
    # classe à splitter 
    ksplit=which.max(varclass)
    ksplit.sauv=c(ksplit.sauv,ksplit)
    
    # splitting, fonctions annexes
    A=FuzzyKmeansSplit(X,mv,gamma,delta,w,ksplit,Pik=probapost[,ksplit],ninit=5,niter=20,K=2)
    
    # mise à jour des centres
    centers[ksplit,]=A$centersNew[1,]
    centers=rbind(centers,A$centersNew[2,])
    nbcluster=nbcluster+1
    rownames(centers)=1:nbcluster
    
    #mise a jour des probapost
    probapost[,ksplit] = A$Pinew[,1]
    probapost=cbind(probapost,A$Pinew[,2])
    
    #mise à jour des poids alpha
    w=weightprobapostbis(X, mv,centers,probapost,gamma,delta)
    weights.sauv[[iter+1]]=w$weights
    
    #mise à jour des variances pondérées par cluster
    varclass=rep(0,nbcluster)
    for (k in 1:nbcluster){
      for (i in 1:nrow(X)){
        for (v in 1:length(mv)){
          dv=(ref[v]+1):ref[v+1]
          varclass[k] = varclass[k] + ( (w$weights[k,v]^gamma) * (probapost[i,k]^delta) * sum( (as.matrix(X)[i,dv] - centers[k,dv])^2))
        }
      }
    }
   
    # calcul du critere general 
    CRIT<-c(CRIT, sum(varclass))
    # 
    iter<-iter+1
  }
return(list(weights=weights.sauv,CRIT=CRIT,withinss=varclass,probapost=probapost,ksplit=ksplit.sauv))  
} 



####################################################
##  Même fonction mais avec les poids alpha sont communs pour les clusters (\alpha_v)
##

splittingProbapost=function(X,mv,gamma=2,delta=2,Kmax,probapost.init){
  #mettre ici un stop si delta <=1 ....
  ref=c(0,cumsum(mv))
  CRIT<-NULL    #evolution critere general
  probapost=probapost.init
  Kinit=ncol(probapost.init)
  
  # initialisation des centres
  centers<-matrix(0,nrow=Kinit,ncol=ncol(X))
  for (k in 1:Kinit){
    centers[k,]<- apply(matrix(rep(probapost[,k]^delta,sum(mv)),nrow=nrow(X)) * X,2,sum) / sum(probapost[,k]^delta)
  }
  rownames(centers)<-1:Kinit
  
  # calcul des poids initiaux
  weights.sauv<-NULL
  w=weightprobapost(X, mv,centers,probapost,gamma,delta)
  weights.sauv<-w$weights
  
  #calcul des variances pondérées par classe   -> optimiser ces boucles for!!
  varclass=rep(0,Kinit)
  for (k in 1:Kinit){
    for (i in 1:nrow(X)){
      for (v in 1:length(mv)){
        dv=(ref[v]+1):ref[v+1]
        varclass[k] = varclass[k] + ( (w$weights[v]^gamma) * (probapost[i,k]^delta) * sum( (as.matrix(X)[i,dv] - centers[k,dv])^2))
      }
    }
  }
  # calcul du critere general 
  CRIT<-c(CRIT, sum(varclass))
  
  ###   début de la boucle de splitting
  ksplit.sauv<-NULL
  iter=1
  nbcluster=Kinit
  
  while(nbcluster < Kmax ){
    # classe à splitter 
    ksplit=which.max(varclass)
    ksplit.sauv=c(ksplit.sauv,ksplit)
    
    # splitting, fonctions annexes
    A=FuzzyKmeansSplit(X,mv,gamma,delta,w,ksplit,Pik=probapost[,ksplit],ninit=20,niter=30,K=2)
    
    # mise à jour des centres
    centers[ksplit,]=A$centersNew[1,]
    centers=rbind(centers,A$centersNew[2,])
    nbcluster=nbcluster+1
    rownames(centers)=1:nbcluster
    
    #mise a jour des probapost
    probapost[,ksplit] = A$Pinew[,1]
    probapost=cbind(probapost,A$Pinew[,2])
    
    #mise à jour des poids alpha
    w=weightprobapost(X, mv,centers,probapost,gamma,delta)
    weights.sauv=cbind(weights.sauv,w$weights)
    
    #mise à jour des variances pondérées par cluster
    varclass=rep(0,nbcluster)
    for (k in 1:nbcluster){
      for (i in 1:nrow(X)){
        for (v in 1:length(mv)){
          dv=(ref[v]+1):ref[v+1]
          varclass[k] = varclass[k] + ( (w$weights[v]^gamma) * (probapost[i,k]^delta) * sum( (as.matrix(X)[i,dv] - centers[k,dv])^2))
        }
      }
    }
    
    # calcul du critere general 
    CRIT<-c(CRIT, sum(varclass))
    # 
    iter<-iter+1
  }
  return(list(weights=weights.sauv,CRIT=CRIT,withinss=varclass,probapost=probapost,ksplit=ksplit.sauv))  
} 









#################################################
# fonction calcul annexe des poids alpha par cluster (\alpha_{k,v})
weightprobapostbis<-function(X, mv,centers,probapost,gamma,delta){
  ref=c(0,cumsum(mv))
  weights=matrix(0,nrow=nrow(centers),ncol=length(mv))
  weightsmv=matrix(0,nrow=nrow(centers),ncol=sum(mv))
  for (k in 1:nrow(centers)){
    aux=NULL
    for (v in 1:length(mv)){
      dv=(ref[v]+1):ref[v+1]
      aux=c(aux, sum((probapost[,k]^delta) * rowSums((sweep(X[,dv,drop=FALSE],2,centers[k,dv],FUN="-"))^2)) )
    }
    if (gamma >1){
      aux=aux^(1/(1-gamma))
      aux=aux/sum(aux)
    }
    if (gamma==1){
      aux=aux^(-1)
      aux=aux/sum(aux)
    }
    aux1<-NULL
    for (j in 1:length(mv))
      aux1<-c(aux1,rep(aux[j],mv[j]))
    weights[k,]=aux
    weightsmv[k,]=aux1
  }  
  return(list(weights=weights,weightsmv=weightsmv))
} 


##################################
# fonction calcul annexe des poids alpha communs aux clusters (alpha_v)
weightprobapost<-function(X, mv, centers, probapost, gamma, delta){
  ref=c(0,cumsum(mv))
  aux=matrix(0,nrow=ncol(probapost),ncol=length(mv))
  for (k in 1:ncol(probapost)){   # for each cluster k
    for (v in 1:length(mv)){
      dv=(ref[v]+1):ref[v+1]
      aux[k,v] =sum( (probapost[,k]^delta) * rowSums((sweep(X[,dv,drop=FALSE],2,centers[k,dv],FUN="-"))^2)) 
    }
  }
  aux1=apply(aux,2,sum)
  
  if (gamma >1){
    aux1=aux1^(1/(1-gamma))
    aux1=aux1/sum(aux1)
  }
  if (gamma==1){
    v0=which.min(aux1)
    aux1[v0]=1
    aux1[-v0]=0
  }
  
  weightsmv<-NULL
  for (j in 1:length(mv))
    weightsmv<-c(weightsmv,rep(aux1[j],mv[j]))
  return(list(weights=aux1,weightsmv=weightsmv))
}  








#####Fuzzy K-means avec vecteur de proba
###  j'ai du mettre plein d'exceptions à cause du cas delta=1 qui mene à des CRIT=NaN
###  la programmation est à améliorer!!!
FuzzyKmeansSplit=function(X,mv,gamma,delta,w,ksplit,Pik,ninit=20,niter=20,K=2){
  #pour stocker les réponses
  Pinew = matrix(0,nrow=nrow(X),ncol=K)
  centersNew=matrix(0,nrow=K,ncol=ncol(X))
  A1=FuzzyKmeansaux(X,w,gamma,delta, ksplit, Pik, K,niter)
  for (z in 1:ninit) {   #boucle sur le nombre d'essai du fuzzy Kmeans
    A2=FuzzyKmeansaux(X,w,gamma,delta, ksplit, Pik, K,niter)
    if(is.na(A1$CRIT)==FALSE){
      if (is.na(A2$CRIT)==FALSE){
        if (A2$CRIT<A1$CRIT){
        A1=A2
        }
      }  
    }else{ 
      if(is.na(A2$CRIT)==FALSE){
        A1=A2
      }
    }
  }  
  Pinew=A1$Pinew
  centersNew=A1$centersNew
return(result=list(Pinew= Pinew, centersNew= centersNew))
}


##  Algo de Fuzzy Kmeans sur les données
FuzzyKmeansaux<-function(X,w,gamma,delta, ksplit, Pik, K,niter){
  
  # Z_i^{(v)} = alpha_v^{gamma/2} X_i^{(v)}
  if (is.vector(w$weights)==TRUE){       # même poids par classe (\alpha_v)
    Z=data.frame(  matrix(rep(w$weightsmv^(gamma/2),nrow(X)),nrow=nrow(X),byrow=T) * X  )
  }else{       #poids(\alpha_{k,v})
    Z=data.frame(  matrix(rep(w$weightsmv[ksplit,]^(gamma/2),nrow(X)),nrow=nrow(X),byrow=T) * X  )  
  }
  
  #init de stockage des Pinew et centersNew
  Pinew=matrix(0,nrow=nrow(X),ncol=K)
  centersNew=matrix(0,nrow=K,ncol=ncol(X))
  eta=matrix(0,nrow=K,ncol=ncol(X))
  m=matrix(0,nrow=nrow(X),ncol=K)   # \|Z_i - \eta_k\|^2
  CRIT<-NULL
  
  #initialisation des centres
  centersNew=kmeans(X,K)$centers
  
  if (delta>1){
    for(l in 1:niter){
      # Mise à jour des poids
      for (k in 1:K){  
        if (is.vector(w$weights)==TRUE){       # même poids par classe (\alpha_v)
          eta[k,] = centersNew[k,] * w$weightsmv^(gamma/2)
        }else{
          eta[k,] = centersNew[k,] * w$weightsmv[ksplit,]^(gamma/2)  
        }
        m[,k]=rowSums((matrix(rep(eta[k,],nrow(Z)),byrow=T,nrow=nrow(Z))-Z)^2)
      }
      Pinew = (m^(1/(1-delta))) / matrix(rep(apply(m^(1/(1-delta)),1,sum),K),ncol=K)
      Pinew =  matrix(rep(Pik,K),ncol=K) * Pinew
      
      # Mise à jour des centres
      for (k in 1:K)
        centersNew[k,] = apply(matrix(rep((Pinew[,k]^delta),ncol(X)),nrow=nrow(X)) * X,2,sum) / sum(Pinew[,k]^delta)
    } # fin de la boucle ninter
    
    # calcul du critere
    for (k in 1:K){ 
      if (is.vector(w$weights)==TRUE){       # même poids par classe (\alpha_v)
        eta[k,] = centersNew[k,] * w$weightsmv^(gamma/2)
      }else{
        eta[k,] = centersNew[k,] * w$weightsmv[ksplit,]^(gamma/2)  
      }
      m[,k]=rowSums((matrix(rep(eta[k,],nrow(Z)),byrow=T,nrow=nrow(Z))-Z)^2)
    }
    CRIT=sum((Pinew^delta) * m)
    
  }else{cat("Error delta value")}
  
  #if(delta==1){
  #  for (l in 1:niter){
  #    for (k in 1:K){  
  #      eta[k,] = centersNew[k,] * w$weightsmv[ksplit,]^(gamma/2)
  #      m[,k]=rowSums((matrix(rep(eta[k,],nrow(Z)),byrow=T,nrow=nrow(Z))-Z)^2)
  #    }
  #    k0=apply(m,1,which.max)
  #    Pinew=matrix(0,nrow=nrow(X),ncol=K)
  #    for (i in 1:nrow(X))
  #      Pinew[i,k0[i]] = Pik[i]
  #      #Pinew = matrix(rep(Pik,K),ncol=K)  *  ((m-matrix(rep(apply(m,1,max),K),ncol=K))<0) 
  #    for (k in 1:K)
  #      centersNew[k,] = apply(matrix(rep((Pinew[,k]^delta),ncol(X)),nrow=nrow(X)) * X,2,sum) / sum(Pinew[,k]^delta)
  #  }  # fin de la boucle iter
  #  for (k in 1:K){  
  #    eta[k,] = centersNew[k,] * w$weightsmv[ksplit,]^(gamma/2)
  #    m[,k]=rowSums((matrix(rep(eta[k,],nrow(Z)),byrow=T,nrow=nrow(Z))-Z)^2)
  #  }
  #  CRIT=sum((Pinew^delta) * m) 
  #  
  #}else{cat("error delta")}
  
  return(list(Pinew=Pinew,centersNew=centersNew,CRIT=CRIT))  
}

# 
# #######################
# ##   exploitation de l'historique du splitting 
# ##  la matrice Msplit explique sur chaque ligne la répartition des classes 
# ##  (1ere ligne le nb classes totales 1 à Kmax   jusqu'à sur la dernière ligne la répartition en 1:Kinit classes)
# 
# MatrixKsplit<-function(ksplit,Kmax){
#   Kref=Kmax
#   Msplit=matrix(rep(1:Kref,2),nrow=2,byrow=T)
#   iter=length(ksplit)
#   while(iter>0){
#   I=c(Kref,which(Msplit[nrow(Msplit)-1,]==Kref))
#   Msplit[nrow(Msplit),I] = ksplit[iter]
#   iter=iter-1
#   Kref=Kref-1
#   Msplit=rbind(Msplit,Msplit[nrow(Msplit),])
#   }
#   return(Msplit[-nrow(Msplit),])
# }
# 
# 
# ########################
# ##   Extraction des probapost dans une étape de l'arbre du splitting
# ##
# ##   Entrée : la matrice finale des probapost, le ksplit et Kref=le nb de classes que l'on souhaite
# 
# ProbapostKref<-function(probapost,ksplit,Kref){
#   Kmax=ncol(probapost)
#   Msplit=MatrixKsplit(ksplit,Kmax)
#   I=which(apply(Msplit,1,max)==Kref)    # la ligne de Msplit qui contient la répartition des Kmax classes en Kref classes
#   probapostextract=matrix(0,nrow=nrow(probapost),ncol=Kref)
#   for (k in 1:Kref)
#     probapostextract[,k]=apply(probapost[,which(Msplit[I,]==k),drop=F],1,sum)
#   
#   return(probapostextract)
# }
# 
# ########################
# ##  Extraction classification par MAP complete ou "smooth"
# ##
# ##  Entrée : matrice de probapost, seuil (par défaut 0.8)
# 
# MAPsmooth<-function(probapost,seuil=0.8){
#   if (seuil<1 & 0<=seuil){
#   label=rep(0,nrow(probapost))
#   label=apply(probapost,1,which.max) * (apply(probapost,1,max)>seuil) 
#   return(label)
#   }
#   else{cat("Error for threshold value")}
# }
# 
# 
