################################################################################
# Loading the packages
################################################################################
library("ape")#
library("apTreeshape")#
library("RColorBrewer")#



################################################################################
# Functions
################################################################################

#### [Create Phylogenies]
phylogenies=function(tree_number,isolation_time,migration,tcorte,tmax,riq_max,filtro){
  

  n<-tree_number;
  l<-isolation_time;
  m<-migration;
  p=200
  dados=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_MS",n,"_.txt"),header=F)#sp, sp recente, sp ancestral, tempo de surgimento, tempo de extincao/sobrevivencia
  obj=read.table(paste0("Pop",p,"_Isol",l,"_Mig",m,"_Temporal",n,"_.txt"),header=TRUE)

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
    if(tn>filtro){
      
      ance=dados[n,2]
      somaAn=sum(M[,ance])
      somaN=sum(M[,n])
      no=no+1
      if(somaAn==0&somaN==0){  
        M_nos[no,1]=paste0("n",no)
        tn=dados[n,4]-dados[n,3]
        tance=dados[ance,4]-dados[n,3]
        M_nos[no,2]=paste0("(",n,":",tn,",",ance,":",tance,")")
        M[no,n]=1
        M[no,ance]=1
        M_nos[no,3]=dados[n,3]
        M_nos$nos[no]=list(c(ance,n))

      }
      if(somaAn>0&somaN>0){
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
      if(somaAn==0&somaN>0){  
        for(nn in 1:(no-1)){
          if(M[nn,n]>0){
            M[no,]=M[no,]+M[nn,]
            noN=M_nos[nn,2]
            nunoN=nn
          }
        }
        M[no,ance]=1
        M_nos[no,1]=paste0("n",no)
        tn=as.numeric(M_nos[nunoN,3])-dados[n,3] 
        tance=dados[ance,4]-dados[n,3]
        M_nos[no,2]=paste0("(",noN,":",tn,",",ance,":",tance,")")
        M_nos[no,3]=dados[n,3]
        M_nos$nos[no]=list(c(ance,M_nos$nos[nunoN][[1]]))
      }
      if(somaAn>0&somaN==0){
        for(nn in 1:(no-1)){
          if(M[nn,ance]>0){ 
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

  #### Achando os edges (t1)
  Edges=t1$edge
  tips=as.numeric(t1$tip.label)
  Tips_Nos=c(as.list(tips),M_nos$nos)
  M_edges=as.data.frame(Edges)
  colnames(M_edges)=c("older","newer") 
  M_edges$cor=0
  for (n in 1:nrow(Edges)){
    esp=Tips_Nos[[n]]
    ed=which.edge(t1,as.character(esp))
    if(length(ed)==1){
      M_edges$cor[ed]=esp 
    }else{
      M_sub=M_edges[ed,]
      ed_2=unique(M_sub$older)
      esp_2=sort(M_sub$cor)[1]
      M_edges$cor[M_edges$newer%in%ed_2&M_edges$cor==0]=esp_2
    }
  }
  
  ### Complete Plot
  first_event = dados[2,3]/2000
  par(mai = c(0.2, 0.05, 0.3, 0.1), oma = c(5, 5, 0.5, 0.5),las=0)
  par(fig=c(first_event,0.97,0.50,1))
  plot(t1,edge.width=2, type = "p", edge.color=cores[M_edges$cor], show.tip.label = F)
  mtext("Complete \n Phylogeny",side=2.2,cex=1.2,line=1,adj=0.5)

  
  #### Achando os edges (t2)
  Edges=t2$edge
  tips=as.numeric(t2$tip.label)
  Tips_Nos=c(as.list(tips),M_nos$nos)
  M_edges=as.data.frame(Edges)
  colnames(M_edges)=c("older","newer") 
  M_edges$cor=0
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

  for (n in nrow(Edges):1){
    if(M_edges$cor[n]==0){
      
      M_sub=M_edges[M_edges$older==M_edges$newer[n],]
      M_sub<-M_sub[order(-M_sub$newer),]
      ed_2=unique(M_sub$older)
      esp_2=M_sub$cor[1]
      M_edges$cor[n]=esp_2
    }
  }
  
  ### Extant Plot
  Anc<-M_edges[order(M_edges$cor),]
  ancestral<-Anc$cor[1]
  first_event = dados[ancestral,3]/2000
  par(fig=c(first_event, 0.97,0.25,0.50), new=T)
  plot(t2,edge.width=2, type = "p", edge.color=cores[M_edges$cor], show.tip.label = F)
  mtext("Extant \nPhylogeny",side=2.2,cex=1.2,line=1,adj=0.5)
  
  
  ### Patches plot
  par(fig=c(0, 0.97,0,0.25),new=T)

  objm=obj[obj$tem<=tmax,]
  objm$nesp=objm$esp
  objm$abund1=(objm$abund1/objm$pop)*0.8/0.5
  objm$abund2=(objm$abund2/objm$pop)*0.8/0.5
  tempos=unique(objm$tem)
  espt=max(unique(objm$nesp))
  times=ceiling(espt/ncor)
  cor=rep(cores, times)
  
  plot(1, type="n", xlab=" ", ylab=" ", xlim=c(tcorte,tmax)
       , ylim=c(0.50, 2.5),yaxt = "n",xaxt = "n")
  axis(2, at=c(1,2), labels=c("P1","P2"),cex.axis=1.0)
  axis(1, at=seq(0,tmax,floor(tmax/4)), labels = F)
  migsv=length(unique(objm$mig))
  mig_min=min(unique(objm$mig))
  
  ti=0
  for(t in tempos){
    i1=objm[objm$tem==t&objm$abund1>0,c("nesp","abund1")]
    i1=i1[order(i1$nesp),]
    mig=unique(objm$mig[objm$tem==t])
    if(migsv>1){
      if(mig==mig_min){
        points(t,1.5, cex=0.5, pch=15)
      }
    }
    y1i=0.6
    for (e in 1:nrow(i1)){
      polygon(x = c(ti,ti,t,t),
              y = c(y1i+i1$abund1[e],y1i,y1i,y1i+i1$abund1[e]),
              density = NA,col = cor[i1$nesp[e]],border = NA) 
      
      y1i=y1i+i1$abund1[e]
    }
    i2=objm[objm$tem==t&objm$abund2>0,c("nesp","abund2")]
    i2=i2[order(i2$nesp),]
    y2i=1.6
    for (e in 1:nrow(i2)){
      polygon(x = c(ti,ti,t,t),
              y = c(y2i+i2$abund2[e],y2i,y2i,y2i+i2$abund2[e]),
              density = NA,col = cor[i2$nesp[e]],border = NA) 
      
      y2i=y2i+i2$abund2[e]
    } 
    ti=t
  }
  

  axis(1, at=seq(0,tmax,floor(tmax/4)),cex.axis=1.5,las=0
       ,labels=as.character(seq(0,2000,500)))
  mtext("Iterations",side=1,cex=1.5,line=4,adj=0.5) 
  
}



# Data repository   #### 
setwd("~/Data/Phylogenies")


      #### ############################################## TO PLOT
      set.seed(018) # just to make it reproducible
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      cores = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      ncor=length(cores)
      cores<-colorRampPalette(cores)(225)
      cores<-unique(cores)
      
      ## Provide data parameters, for example:
      tree_number=1
      isolation_time=80
      migration=0.002
      dev.off()
      
      # Print phylogenies
      phylogenies(tree_number
                  ,isolation_time
                  ,migration
                  ,0
                  ,2000
                  ,20
                  ,0  # filtro geracional (under construction)
                  ) 

################################################################################

### ############################################################################

      # Heavly Modified in: 04/05/23 ####

