################################################################################
#loading the packages
################################################################################
library("ape")
library("vegan")
#library("phytools")
library("ggplot2")
library("cowplot") 
library("apTreeshape")
library("gridExtra")
library("viridis")
library("RColorBrewer")
library("ggcorrplot")
library("corrplot")
library("reshape")
library("tidyr")
library("tidyverse")


################################################################################
#analises de estrutura de filogenia
################################################################################

#### [Gamma-statistics]
gamma_steps<-function(Nspp,alpha0){
  gamma_values=data.frame(0,0)
  colnames(gamma_values)=c("gamma0","gammap0")
  
  tk=rep(0,Nspp)
  tkp=rep(0,Nspp)
  tk = 0.0
  tkp = 0.0
  tk[2] = 2.0^(1.0-alpha0)
  tkp[2] = -log(2.0)*tk[2]
  for(i in 3:Nspp){
    ai = as.double(i)
    aux = ai^(1.0-alpha0)
    tk[i] = tk[i-1] + aux
    tkp[i] = tkp[i-1] - log(ai)*aux
  }
  
  deno = tk[Nspp]/sqrt(as.double(Nspp-2)*12.0)
  
  aux1 = 0.0
  aux2 = 0.0
  for(i in 2:(Nspp-1)){
    aux1 = aux1 + tk[i]
    aux2 = aux2 + tkp[i]
  }
  
  gamma_values$gamma0 = (aux1/as.double(Nspp-2) - 0.5*tk[Nspp])/deno
  gamma_values$gammap0 = (aux2 - aux1*tkp[Nspp]/tk[Nspp])/(deno*as.double(Nspp-2))
  
  return(gamma_values)
}
#### [Alpha-value]
alpha_value <- function(Nspp,gm){
  error = 0.001
  alpha0 = 1.0
  deltaalpha = 1.0
  
  while(abs(deltaalpha) > error){
    gm_calc=gamma_steps(Nspp,alpha0)
    deltaalpha = (gm - gm_calc$gamma0)/gm_calc$gammap0
    alpha0 = alpha0 + deltaalpha
  }
  alpha = alpha0
  return(alpha)
}
#### [Phylogeny-age]
phylo_age<-function(dados,obj,tcorte,tmax,riq_max,filtro){
  #filtro<-1
  #dados=read.table("Pop200_Isol100_Mig0.01_MS22_.txt")#sp, sp recente, sp ancestral, tempo de surgimento, tempo de extincao/sobrevivencia
  dim=nrow(dados)-1
  M=matrix(0,dim,dim+1)
  #M_nos=data.frame(No=as.character(rep("",dim)), Tree=as.character(rep("",dim)), Time=as.numeric(rep(-1,dim)), nos=as.numeric(rep(-1,dim)))
  M_nos=data.frame(No=as.character(), Tree=as.character(), Time=as.numeric(), nos=as.numeric())
  M_nos[1:(dim),]=NA
  M_nos$No="gata"
  M_nos$Tree="branca"
  M_nos$Time=-1
  M_nos$nos=-1
  no=0
  #n=56
  for(n in nrow(dados):2){
    
    tn=dados[n,4]-dados[n,3]
    if(tn>filtro){
      
      ance=dados[n,2]#ancestral cluna 2
      #verica se o ancestral ja esta em outro no
      somaAn=sum(M[,ance])
      somaN=sum(M[,n])
      no=no+1
      if(somaAn==0&somaN==0){  #espécie já foi usada?
        #M_nos[no,] = c("a","b",-1,-1)
        M_nos[no,1]=paste0("n",no)
        #M_nos$No[no]=paste0("n",no)
        tn=dados[n,4]-dados[n,3]# tempo da esp?cie
        tance=dados[ance,4]-dados[n,3]#tance tempo do ancestral
        M_nos[no,2]=paste0("(",n,":",tn,",",ance,":",tance,")")
        M[no,n]=1
        M[no,ance]=1
        M_nos[no,3]=dados[n,3]#tempo
        M_nos$nos[no]=list(c(ance,n))
        #if(tance>tn){M_nos[no,3]=tance}
      }
      if(somaAn>0&somaN>0){ #no que ja foi usado 
        for(nn in 1:(no-1)){
          if(M[nn,n]>0){
            M[no,]=M[no,]+M[nn,]
            noN=M_nos[nn,2]
            nunoN=nn
          }
          
          if(M[nn,ance]>0){
            M[no,]=M[no,]+M[nn,]
            noAnce=M_nos[nn,2]
            nunoAnce=nn
          }
        }
        M_nos[no,1]=paste0("n",no)
        tn=as.numeric(M_nos[nunoN,3])-dados[n,3]
        tance=as.numeric(M_nos[nunoAnce,3])-dados[n,3]
        M_nos[no,2]=paste0("(",noN,":",tn,",",noAnce,":",tance,")")
        M_nos[no,3]=dados[n,3]
        M_nos$nos[no]=list(c(M_nos$nos[nunoAnce][[1]],M_nos$nos[nunoN][[1]]))
      }
      if(somaAn==0&somaN>0){  #será que este if faz sentido?
        for(nn in 1:(no-1)){
          if(M[nn,n]>0){  #resgata info do n
            M[no,]=M[no,]+M[nn,]
            noN=M_nos[nn,2]
            nunoN=nn
          }
        }
        M[no,ance]=1
        M_nos[no,1]=paste0("n",no)
        tn=as.numeric(M_nos[nunoN,3])-dados[n,3] #o numero do nó ancetral -tempo que o no anterior surgiu
        tance=dados[ance,4]-dados[n,3]
        M_nos[no,2]=paste0("(",noN,":",tn,",",ance,":",tance,")")
        M_nos[no,3]=dados[n,3]
        M_nos$nos[no]=list(c(ance,M_nos$nos[nunoN][[1]]))
      }
      if(somaAn>0&somaN==0){
        for(nn in 1:(no-1)){
          if(M[nn,ance]>0){  #resgata info do n
            M[no,]=M[no,]+M[nn,]
            noAnce=M_nos[nn,2]
            nunoAnce=nn
          }
        }
        M[no,n]=1
        M_nos[no,1]=paste0("n",no)
        tn=dados[n,4]-dados[n,3]
        tance=as.numeric(M_nos[nunoAnce,3])-dados[n,3]
        M_nos[no,2]=paste0("(",n,":",tn,",",noAnce,":",tance,")")
        M_nos[no,3]=dados[n,3]
        M_nos$nos[no]=list(c(M_nos$nos[nunoAnce][[1]],n))
      }
    }
  }      
  
  t1=read.tree(text=paste0(M_nos[no,2],";"))
  t2<-drop.fossil(t1)
  
  #### Achando os edges (t2)
  Edges=t2$edge
  tips=as.numeric(t2$tip.label)
  Tips_Nos=c(as.list(tips),M_nos$nos)
  M_edges=as.data.frame(Edges)
  colnames(M_edges)=c("older","newer") 
  M_edges$cor=0
  #M_edges <-M_edges[order(M_edges$older),]
  #n=6
  for (n in 1:nrow(Edges)){
    esp=Tips_Nos[[n]]
    ed=which.edge(t2,as.character(esp))
    if(length(ed)==1){
      M_edges$cor[ed]=esp 
    }else{
      M_sub=M_edges[ed,]
      ed_2=
        esp_2=sort(M_sub$cor)[1]
      M_edges$cor[M_edges$newer%in%ed_2&M_edges$cor==0]=esp_2
    }
  }
  #n=6
  for (n in nrow(Edges):1){
    if(M_edges$cor[n]==0){
      
      M_sub=M_edges[M_edges$older==M_edges$newer[n],]
      M_sub<-M_sub[order(-M_sub$newer),]
      ed_2=unique(M_sub$older)
      esp_2=M_sub$cor[1]
      M_edges$cor[n]=esp_2
    }
  }
  
  Anc<-M_edges[order(M_edges$cor),]
  ancestral<-Anc$cor[1]
  first_event = 2000-dados[ancestral,3]
  
  
  return(first_event)
}
#### [Dirvesity]
diversity<-function(objT){
  
  
  #t=200;l=100;m=0.2;n=20;p=200
  #objT=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_Temporal",n,"_.txt"),header=TRUE)
  #Div=data.frame(tem=as.numeric(), N1=as.numeric(), N2=as.numeric(), NT=as.numeric(), 
  #               Beta=as.numeric(), D0=as.numeric(), D1=as.numeric(), D2=as.numeric())
  Div <- as.data.frame(matrix(0, ncol = 10, nrow = 50))
  colnames(Div)=c("tem","N1","N2","NT","Beta","D0","D1","D2","H","|N1-N2|")
  Div[1:(200),]=NA
  for(t in 0:200){
    
    
    
    
    Beta=objT[objT$tem==t*10,]
    Beta$esp1[Beta$abund1!=0]=1;
    #Beta$esp1[Beta$esp1==NA]=0
    Beta$esp2[Beta$abund2!=0]=1;
    #Beta$esp2[Beta$esp2==NA]=0
    objT$D=(objT$abund1+objT$abund2)/objT$pop
    
    Div[t,1]=t*10
    Div[t,2]=sum(Beta$esp1,na.rm=T)
    Div[t,3]=sum(Beta$esp2,na.rm=T)
    Div[t,4]=length(Beta$tem)
    Div[t,5]=(2*Div[t,4] - Div[t,2] -Div[t,3])/Div[t,4]
    num<- (Div[t,2] -Div[t,3])*(Div[t,2] -Div[t,3])
    den<- (Div[t,2] +Div[t,3])
    #Div[t,10]=abs(Div[t,2] -Div[t,3])
    Div[t,10]=sqrt(num)/den
    
    
    D0=sum((Beta$D)^0,na.rm=T)
    Div[t,6]=D0
    D1=1#exp(-sum((Beta$D*log(Beta$D)),na.rm=T))
    Div[t,7]=D1
    D2=sum((Beta$D)^2,na.rm=T)
    D2=1/D2
    Div[t,8]=D2
    Div[t,9]=1 - 1/D2
   
    
    
    
    
  }
  return(Div)
}

#### [Supersymmetric Trees] (necessary to Gamma and Alpha)
timetree <- function(tree)
  try(chronoMPL(multi2di(tree)), silent=TRUE)

#### [Verify]
verify<-function(dados,filtro){
  
  cont=0
  for(i in nrow(dados):2){
    #print(i)
    tn=dados[i,4]-dados[i,3]
    if(tn>filtro){
      cont<-cont+1
    }
  }
  
  v=cont
  return(v)
}
#### [Criar tree]
createTree<-function(dados,obj,tcorte,tmax,riq_max,filtro){
  #dados=read.table("Isol100_Mig0.16_Temporal1_.txt",header=TRUE)#sp, sp recente, sp ancestral, tempo de surgimento, tempo de extincao/sobrevivencia
  dim=nrow(dados)-1
  M=matrix(0,dim,dim+1)
  M_nos=data.frame(No=as.character(), Tree=as.character(), Time=as.numeric(), nos=as.numeric())
  M_nos[1:(dim),]=NA
  M_nos$No="gata"
  M_nos$Tree="branca"
  M_nos$Time=-1
  M_nos$nos=-1
  no=0
  
  for(n in nrow(dados):2){
    
    
    tn=dados[n,4]-dados[n,3]
    #print(tn)
    #print(n)
    #if((tn>filtro) || (dados[n,4]=2001)){
    if(tn>filtro){
      
      ance=dados[n,2]#ancestral cluna 2
      #verica se o ancestral ja esta em outro no
      somaAn=sum(M[,ance])
      somaN=sum(M[,n])
      no=no+1
      if(somaAn==0&somaN==0){  #no novo
        M_nos[no,1]=paste0("n",no)
        tn=dados[n,4]-dados[n,3]# tempo da esp?cie
        tance=dados[ance,4]-dados[n,3]#tance tempo do ancestral
        M_nos[no,2]=paste0("(",n,":",tn,",",ance,":",tance,")")
        M[no,n]=1
        M[no,ance]=1
        M_nos[no,3]=dados[n,3]#tempo 
        #if(tance>tn){M_nos[no,3]=tance}
      }
      if(somaAn>0&somaN>0){ #no que ja foi usado 
        for(nn in 1:(no-1)){
          if(M[nn,n]>0){
            M[no,]=M[no,]+M[nn,]
            noN=M_nos[nn,2]
            nunoN=nn
          }
          
          if(M[nn,ance]>0){
            M[no,]=M[no,]+M[nn,]
            noAnce=M_nos[nn,2]
            nunoAnce=nn
          }
        }
        M_nos[no,1]=paste0("n",no)
        tn=as.numeric(M_nos[nunoN,3])-dados[n,3]
        tance=as.numeric(M_nos[nunoAnce,3])-dados[n,3]
        M_nos[no,2]=paste0("(",noN,":",tn,",",noAnce,":",tance,")")
        M_nos[no,3]=dados[n,3]
      }
      if(somaAn==0&somaN>0){  #será que este if faz sentido?
        for(nn in 1:(no-1)){
          if(M[nn,n]>0){  #resgata info do n
            M[no,]=M[no,]+M[nn,]
            noN=M_nos[nn,2]
            nunoN=nn
          }
        }
        M[no,ance]=1
        M_nos[no,1]=paste0("n",no)
        tn=as.numeric(M_nos[nunoN,3])-dados[n,3] #o numero do nó ancetral -tempo que o no anterior surgiu
        tance=dados[ance,4]-dados[n,3]
        M_nos[no,2]=paste0("(",noN,":",tn,",",ance,":",tance,")")
        M_nos[no,3]=dados[n,3]
        
      }
      if(somaAn>0&somaN==0){
        for(nn in 1:(no-1)){
          if(M[nn,ance]>0){  #resgata info do n
            M[no,]=M[no,]+M[nn,]
            noAnce=M_nos[nn,2]
            nunoAnce=nn
          }
        }
        M[no,n]=1
        M_nos[no,1]=paste0("n",no)
        tn=dados[n,4]-dados[n,3]
        tance=as.numeric(M_nos[nunoAnce,3])-dados[n,3]
        M_nos[no,2]=paste0("(",n,":",tn,",",noAnce,":",tance,")")
        M_nos[no,3]=dados[n,3]
        
        
      }
    }
  }       
  t=read.tree(text=paste0(M_nos[no,2],";"))
  #suricat
  
  #### Achando os edges
  Edges=t$edge
  tips=as.numeric(t$tip.label)
  Tips_Nos=c(as.list(tips),M_nos$nos)
  M_edges=as.data.frame(Edges)
  colnames(M_edges)=c("older","newer") 
  M_edges$cor=0
  for (n in 1:nrow(Edges)){
    esp=Tips_Nos[[n]]
    ed=which.edge(t,as.character(esp))
    if(length(ed)==1){
      M_edges$cor[ed]=esp 
    }else{
      M_sub=M_edges[ed,]
      ed_2=unique(M_sub$older)
      esp_2=sort(M_sub$cor)[1]
      M_edges$cor[M_edges$newer%in%ed_2&M_edges$cor==0]=esp_2
    }
  }
  
  # PLOT Tree
  #first_event = dados[2,3]/2000
  #par(fig=c(first_event, 1,0.50,1))
  #plot(t,edge.width=2, type = "p", edge.color=cores[M_edges$cor], show.tip.label = F)
  
  
  
  #### SPECIES SIZE OF EACH ISLAND
  #par(fig=c(0, 0.97,0,0.35),new=T)
  #objm=obj[obj$tem<=tmax,]
  #objm$nesp=objm$esp
  #objm$abund1=(objm$abund1/objm$pop)*0.8/0.5
  #objm$abund2=(objm$abund2/objm$pop)*0.8/0.5
  #tempos=unique(objm$tem)
  #espt=max(unique(objm$nesp))
  #times=ceiling(espt/ncor)
  #cor=rep(cores, times)
  
  #plot(1, type="n", xlab=" ", ylab=" ", xlim=c(tcorte,tmax), ylim=c(0.50, 2.5),yaxt = "n",xaxt = "n")
  #axis(2, at=c(1,2), labels=c("island 1","island 2"),cex.axis=1.5)
  #axis(1, at=seq(0,tmax,floor(tmax/4)), labels = F)
  #migsv=length(unique(objm$mig))
  #mig_min=min(unique(objm$mig))
  
  
  
  return(t)
  
  
}
#### [PLOTs] as in (Debora et. al, 2022)


###################################### [J function] as in (Lemant et. al, 2022)
move_up <- function(edges, identity) {
  if(!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) stop("Invalid identity.")
  parent <- edges[which(edges$Identity == identity), "Parent"]
  if(length(parent) == 0) return(identity) # if identity is the root then don't move
  if(is.factor(parent)) parent <- levels(parent)[parent]
  return(parent)
}
find_root <- function(edges) {
  start <- edges$Parent[1] # reasonable guess
  if(is.factor(start)) start <- levels(start)[start]
  repeat {
    if(move_up(edges, start) == start) break
    start <- move_up(edges, start)
  }
  return(start)
}
get_subtree_sizes <- function(tree,i=NULL,Adj=NULL,Col=NULL,Cumul=NULL,is_leaf=NULL){
  n<-length(tree$Identity)
  has_pops <- FALSE
  if("Population" %in% colnames(tree)) has_pops <- TRUE
  if(is.null(Adj)) Adj <- get_Adj(tree)
  if(is.null(i)) i <- which(tree$Identity == find_root(tree[,1:2]))
  if(is.null(Col)) {
    Col <- rep("w",n)
    names(Col) <- unique(tree$Identity)
  }
  if(is.null(Cumul)) {
    Cumul <- rep(NA,n)
    names(Cumul) <- unique(tree$Identity)
  }
  if(is.null(is_leaf)) {
    is_leaf <- rep(FALSE, n)
    names(is_leaf) <- unique(tree$Identity)
  }
  if(is.null(Adj[[i]])) is_leaf[i] <- TRUE
  for (j in Adj[[i]]){
    if (Col[j] == "w"){
      L <- get_subtree_sizes(tree,j,Adj,Col,Cumul,is_leaf)
      Col<- L$colour
      Cumul <- L$cumulative
      is_leaf <- L$is_leaf
    }
  }
  Col[i] <- "b"
  if(has_pops) {
    Cumul[i] <- tree$Population[i] + sum(Cumul[Adj[[i]]])
  } else {
    Cumul[i] <- ifelse(is_leaf[i] == TRUE, 1, 0) + sum(Cumul[Adj[[i]]])
  }
  return(list("colour"=Col,"cumulative"=Cumul,"is_leaf"=is_leaf))
}
get_Adj <- function(tree) {
  n<-length(tree$Identity)
  Adj <- vector(mode = "list", length = n)
  for (i in 1:n) if(tree$Parent[i] != tree$Identity[i]) {
    p <- which(tree$Identity == tree$Parent[i])
    Adj[[p]] <- append(Adj[[p]], i)
  }
  return(Adj)
}
J1_index <- function(tree, q = 1, nonrootdomfactor = FALSE) {
  
  if(!is.na(tree)[1]) {
    if("phylo" %in% class(tree)) { # convert from phylo object to appropriate data frame
      tree <- tree$edge
      tree <- as.data.frame(tree)
      colnames(tree) <- c("Parent", "Identity")
    }
    tree <- na.omit(tree) # remove any rows containing NA
    if(is.factor(tree$Parent)) tree$Parent <- levels(tree$Parent)[tree$Parent]
    if(is.factor(tree$Identity)) tree$Identity <- levels(tree$Identity)[tree$Identity]
    start <- setdiff(tree$Parent, tree$Identity)
    if(length(start) > 0) { # add row for root node
      if("Population" %in% colnames(tree)) {
        root_row <- data.frame(Parent = start, Identity = start, Population = 0)
        #message("Assigning Population = 0 to the root node")
      }
      else root_row <- data.frame(Parent = start, Identity = start)
      tree <- rbind(root_row, tree)
    }
  }
  
  n<-length(tree$Identity)
  if (n<=1) return(0)
  Adj <- get_Adj(tree) # adjacency list
  subtree_sizes <- get_subtree_sizes(tree, Adj = Adj) # get the list of all subtree sizes
  Cumul <- subtree_sizes$cumulative # subtree sizes, including the root
  eff_int_nodes <- which(!subtree_sizes$is_leaf) # vector of internal nodes
  leaves <- which(subtree_sizes$is_leaf) # vector of leaves
  # if population sizes are missing then assign size 0 to internal nodes, and size 1 to leaves:
  if(!("Population" %in% colnames(tree))) {
    tree$Population <- rep(0, n)
    tree$Population[leaves] <- 1
    #message("Assigning Population = 0 to internal nodes and Population = 1 to leaves")
  }
  J <- 0
  Star <- Cumul - tree$Population # subtree sizes, excluding the root
  for (i in 1:n){ # loop over all nodes
    if (Star[i] > 0){ # if node has at least one child with non-zero size
      K <- 0
      if(length(Adj[[i]])>1){ # otherwise i has only one child and its balance score is 0
        eff_children <- 0 # number of children with non-zero size
        for (j in Adj[[i]]){
          if (Cumul[j]>0){ # otherwise child j has a 0-sized subtree and does not count
            eff_children <- eff_children+1
            # p is the ratio of the child subtree size including the root (root = the child) 
            # to the parent subtree size excluding the root
            p <- Cumul[j]/Star[i]
            # K is the sum of the node balance scores
            if(q == 1) {
              K <- K + -p*log(p)
            } else {
              K <- K + p^q
            }
          }
        }
        # non-root dominance factor:
        if(nonrootdomfactor) {
          h_factor <- Star[i] / Cumul[i]
        } else {
          h_factor <- 1
        }
        # normalize the sum of balance scores, adjust for non-root dominance, 
        # and then add the result to the index
        if(q == 1) {
          J <- J + h_factor * Star[i] * K / log(eff_children)
        } else {
          J <- J + h_factor * Star[i] * (1 - K) * eff_children^(q - 1) / (eff_children^(q - 1) - 1)
        }
      }
    }
  }
  # normalize the index by dividing by the sum of all subtree sizes:
  if (length(eff_int_nodes)>0) J <- J/sum(Star[eff_int_nodes])
  return(as.numeric(J))
}

###################################### Functions to reproduce the labeled figures

################################################################################








# GENERATE YOUR OWN OUTPUT DATA   #### 

# Set the directory containing the phylogenies
setwd("~/Data/Phylogenies")

# Total number of trees for each combination of parametes
num_trees<-50
# Isolation time in level list
level<-c(80,40,20,0)
# Mean migration rate in mig list
mig<-c(0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.08)

# Generating .RData for complete trees and correlation values
complete<-function(num_trees){
  
  mat    <- as.data.frame(matrix(0, ncol = length(mig), nrow = length(level), dimnames=list(level, mig)))
  Mean_comp <- lapply(seq_len(7), function(x) mat)
  Inf_comp  <- lapply(seq_len(7), function(x) mat)
  Sup_comp  <- lapply(seq_len(7), function(x) mat)
  
  output_values <- as.data.frame(matrix(0, ncol = 7, nrow = num_trees))
  correlation_values <- as.data.frame(matrix(0, ncol =9, nrow = num_trees*length(mig)*length(level)))
  Na <- as.data.frame(matrix(0, ncol = 7, nrow = 1))
  colnames(correlation_values)=c("mig","level","J","Age","Alpha","Asym","Beta","Richness","Speciation")
  colnames(output_values)=c("total_spec","rich","Beta","Asym","J","alpha","Age")
  colnames(Na)=c("total_spec","rich","Beta","Asym","J","alpha","Age")
  
  
  filtro<-0
  y<-0    
  c<-0
  p<-200
  for(l in level){ 
    y<-y+1
    x<-0
    for(m in mig){
      x<-x+1
      na<-0
      
      for(n in 1:num_trees){
        c<-c+1
        dados=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_MS",n,"_.txt"),header=F)#sp, sp recente, sp ancestral, tempo de surgimento, tempo de extincao/sobrevivencia
        obj1=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_Temporal",n,"_.txt"),header=TRUE)
        rep=1
        objT=obj1[obj1$rep==rep,]
        
        d<-diversity(objT)
        Beta<-d[d$tem==2000,5]
        Asym<- d[d$tem==2000,10]
        TajD=1
        
        #### Species Counting
        # Total number of speciation events
        total_spec<-nrow(dados)
        # Total number of extinction events
        cont<-0
        for(s in nrow(dados)){
          if(dados[s,4]==2001){
            cont<-cont+1
          }
        }
        total_ext<-nrow(dados)-cont
        # Species Richness at the end of process
        sub=objT[objT$tem==2000,]
        rich<-nrow(sub)
        
        
        #### Phylogeny metrics
        if (nrow(dados)>1){
          v<-verify(dados,filtro)
          if(v>0){
            
            a<-createTree(dados,objT,0,2000,20,filtro)
            Age<-2000-dados[2,3]
            #calculating J_index
            J1<-J1_index(a)
            #number of spp/tips
            Ntips=length(a$tip.label) 
            
            if(Ntips>2){
              
              #supersymmetry transformation
              Phylo<-timetree(a)
              #calculating gamma
              Gamma=gammaStat(Phylo) 
              #calculating alpha
              Alpha=alpha_value(Ntips,Gamma)
              #calculating sackin
              treeShape<-as.treeshape(Phylo)
              sackin<-sackin(treeShape,norm="yule")
            } else{Alpha=NA;Gamma=NA;}
          } else{Alpha=NA;Gamma=NA;}
        }else{Alpha=NA;Gamma=NA;}
        
        
        
        myList <- list(total_spec,rich,Beta,Asym,J1,Alpha,Age)
        for (i in 1:length(myList)) {
          output_values[n,i]<-myList[i]
        }
        for (i in 1:length(myList)) {
          Na[1,i]<-sum(is.na(output_values[,i]))
        }
        myList2 <- list(m,l,J1,Age,Alpha,Asym,Beta,rich,total_spec) 
        #c("mig","level","J","Age","Alpha","Asym","Beta","Richness","Speciation")
        for (i in 1:length(myList2)) {
          correlation_values[c,i]<-myList2[i]
        }
        
        
        
      }   
      
      #### If at least 20% are different from NA: group the data mean values and confidence interval
      # Mean_list <- list(spe,rich,Beta,Asym,J,alp,Age)   
      for (i in 1:length(Mean_comp)) {
        if(Na[1,i]<=40){
          
          q<-quantile(output_values[,i], prob = c(0.05, 0.95), na.rm = TRUE)
          Sup_comp[[i]][y,x]  <-q[2]
          Inf_comp[[i]][y,x]  <-q[1]
          Mean_comp[[i]][y,x] <-mean(output_values[,i], na.rm=TRUE) 
          
        }else{
          
          Sup_comp[[i]][y,x]  <-NA
          Inf_comp[[i]][y,x]  <-NA
          Mean_comp[[i]][y,x] <-NA
          
        }
      }
      
      
    }
  }
  
  # COMPLETE TREE DATA
  saveRDS(Mean_comp, file="Mean_comp.RData")
  saveRDS(Inf_comp, file="Inf_comp.RData")
  saveRDS(Sup_comp, file="Sup_comp.RData")
  
  # CORRELATION
  saveRDS(correlation_values, file="correlation_values.RData")
}
complete(num_trees)

# Generating .RData for extant trees
extant<-function(num_trees){
  
  mat    <- as.data.frame(matrix(0, ncol = length(mig), nrow = length(level), dimnames=list(level, mig)))
  Mean_ext <- lapply(seq_len(7), function(x) mat)
  Inf_ext  <- lapply(seq_len(7), function(x) mat)
  Sup_ext  <- lapply(seq_len(7), function(x) mat)
  
  output_values <- as.data.frame(matrix(0, ncol = 7, nrow = num_trees))
  correlation_values <- as.data.frame(matrix(0, ncol =9, nrow = num_trees*length(mig)*length(level)))
  Na <- as.data.frame(matrix(0, ncol = 7, nrow = 1))
  colnames(correlation_values)=c("mig","level","J","Age","Alpha","Asym","Beta","Richness","Speciation")
  colnames(output_values)=c("total_spec","rich","Beta","Asym","J","alpha","Age")
  colnames(Na)=c("total_spec","rich","Beta","Asym","J","alpha","Age")
  
  
  filtro<-0
  y<-0    
  c<-0
  p<-200
  for(l in level){ 
    y<-y+1
    x<-0
    for(m in mig){
      x<-x+1
      na<-0
      
      for(n in 1:num_trees){
        c<-c+1
        dados=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_MS",n,"_.txt"),header=F)#sp, sp recente, sp ancestral, tempo de surgimento, tempo de extincao/sobrevivencia
        obj1=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_Temporal",n,"_.txt"),header=TRUE)
        rep=1
        objT=obj1[obj1$rep==rep,]
        
        d<-diversity(objT)
        Beta<-d[d$tem==2000,5]
        Asym<- d[d$tem==2000,10]
        TajD=1
        
        #### Species Counting
        # Total number of speciation events
        total_spec<-nrow(dados)
        # Total number of extinction events
        cont<-0
        for(s in nrow(dados)){
          if(dados[s,4]==2001){
            cont<-cont+1
          }
        }
        total_ext<-nrow(dados)-cont
        # Species Richness at the end of process
        sub=objT[objT$tem==2000,]
        rich<-nrow(sub)
        
        
        #### Phylogeny metrics
        if (nrow(dados)>1){
          v<-verify(dados,filtro)
          if(v>0){
            
            a<-createTree(dados,objT,0,2000,20,filtro)
            a<-drop.fossil(a)
            Age<-phylo_age(dados,objT,0,2000,20,filtro)
            #calculating J_index
            J1<-J1_index(a)
            #number of spp/tips
            Ntips=length(a$tip.label) 
            if(rich==1){Age=NA;J1=NA;}
            if(Ntips>2){
              
              #supersymmetry transformation
              Phylo<-timetree(a)
              #calculating gamma
              Gamma=gammaStat(Phylo) 
              #calculating alpha
              Alpha=alpha_value(Ntips,Gamma)
              #calculating sackin
              treeShape<-as.treeshape(Phylo)
              sackin<-sackin(treeShape,norm="yule")
            } else{Alpha=NA;Gamma=NA;}
          } else{Alpha=NA;Gamma=NA;}
        }else{Alpha=NA;Gamma=NA;}
 
        
        
        myList <- list(total_spec,rich,Beta,Asym,J1,Alpha,Age)
        for (i in 1:length(myList)) {
          output_values[n,i]<-myList[i]
        }
        for (i in 1:length(myList)) {
          Na[1,i]<-sum(is.na(output_values[,i]))
        }
      }   
      
      #### If at least 20% are different from NA: group the data mean values and confidence interval
      # Mean_list <- list(spe,rich,Beta,Asym,J,alp,Age)  
      for (i in 1:length(Mean_ext)) {
        if(Na[1,i]<=40){
          
          q<-quantile(output_values[,i], prob = c(0.05, 0.95), na.rm = TRUE)
          Sup_ext[[i]][y,x]  <-q[2]
          Inf_ext[[i]][y,x]  <-q[1]
          Mean_ext[[i]][y,x] <-mean(output_values[,i], na.rm=TRUE) 
          
        }else{
          
          Sup_ext[[i]][y,x]  <-NA
          Inf_ext[[i]][y,x]  <-NA
          Mean_ext[[i]][y,x] <-NA
          
        }
      }
      
      
    }
  }
  
  
  # EXTANT TREE DATA
  saveRDS(Mean_ext, file="Mean_ext.RData")
  saveRDS(Inf_ext, file="Inf_ext.RData")
  saveRDS(Sup_ext, file="Sup_ext.RData")
  
}
extant(num_trees)

# J and Alpha over evolutionary time
# Provide the total number of trees (num_trees), and the corresponding isolation time (l) and migration rate (m)  
over_time(0.08,80,50)
over_time<-function(migration, isolation_time, num_trees){
  
  #### Estimated mean waiting time to first speciation event (iteration = 260)
  time <- seq(260,2000,10)
  output <- as.data.frame(matrix(0, ncol = 2, nrow = 50))
  to_plot <- as.data.frame(matrix(0, ncol = 7, nrow = length(time)))
  colnames(to_plot)=c("time","J","Jinf","Jsup","Alpha","Ainf","Asup")
  
  
  
  l=isolation_time; m=migration; p=200
  for(j in 1:length(time))  {
    to_plot[j,1]<-time[j]
    for(n in 1:num_trees){ 
      dados=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_MS",n,"_.txt"),header=F)#sp, sp recente, sp ancestral, tempo de surgimento, tempo de extincao/sobrevivencia
      objT=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_Temporal",n,"_.txt"),header=TRUE)
      
      dv<-subset(dados,dados$V3<=time[j])
      if(nrow(dv)>2)
      {
        for(i in row(dv)){
          if(dv[i,4]>time[j]){
            dv[i,4]<-time[j]
          }
        }
        
        a<-createTree(dv,objT,0,2000,20,0)
        #a<-drop.fossil(a) #### Uncomment to extant tree
        Ntips=length(a$tip.label)
        if(Ntips>=2){ output[n,1]<-J1_index(a)}else{output[n,1]<-NA}
        
        Ntips=length(a$tip.label)
        if(Ntips>2){
          Phylo<-timetree(a)
          Gamma=gammaStat(Phylo) 
          output[n,2]<-alpha_value(Ntips,Gamma)
        } else{output[n,2]<-NA}
        
      }
    }
    
    #colnames(to_plot)=c("time","J","Jinf","Jsup","Alpha","Ainf","Asup")
    q<-quantile(output[,1], prob = c(0.05, 0.95), na.rm = TRUE)
    to_plot[j,2] <-mean(output[,1], na.rm=TRUE) 
    to_plot[j,3]   <-q[1]
    to_plot[j,4]   <-q[2]
    
    q<-quantile(output[,2], prob = c(0.05, 0.95), na.rm = TRUE)
    to_plot[j,5] <-mean(output[,2], na.rm=TRUE) 
    to_plot[j,6]   <-q[1]
    to_plot[j,7]   <-q[2]
  }
  
  saveRDS(to_plot, file="over_time.RData")
  
}


################################################################################


# Heavly Modified in: 04/05/23 ####


