# Loading packages
library("reshape")
library("ggplot2")
library("gridExtra")
library("corrplot")
library("RColorBrewer")

###################################### Functions to reproduce the labeled figures
fig3<-function(){
  
  pcomp<-plot_function(Comp_mean, Comp_inf, Comp_sup)
  pcomp[[1]]<-pcomp[[1]]+labs( x =" ") + fundo_branco_grid +  ggtitle("a)")
  pcomp[[2]]<-pcomp[[2]]+labs( x =" ") + fundo_branco_grid + ggtitle("b)") 
  pcomp[[3]]<-pcomp[[3]]+labs( x = expression(bar(epsilon))) + fundo_branco_grid + ggtitle("c)") 
  pcomp[[4]]<-pcomp[[4]]+labs( x =expression(bar(epsilon))) + fundo_branco_grid + ggtitle("d)") 
  plots <- list(pcomp[[1]],pcomp[[2]],
                pcomp[[3]],pcomp[[4]])
  
  legend<-pcomp[[1]]+theme(legend.position = "bottom",
                           legend.title = element_text(size=20),
                           legend.text  = element_text(size=20))
  shared_legend<-get_legend(legend)
  
  grobs <- list()
  widths <- list()
  for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  
  maxwidth <- do.call(grid::unit.pmax, widths)
  
  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  dev.off()
  
  
  #setwd("~/Documents")
  #pdf("fig3.pdf",width = 10,height=6.5)
  p <-grid.arrange(arrangeGrob(
    grobs[[1]],grobs[[2]],
    grobs[[3]],grobs[[4]],nrow = 2, ncol = 2),
    shared_legend, nrow=2, heights = c(10, 1))
  #dev.off()
  
}
fig4<-function(){
  
  pext<-plot_function(Ext_mean, Ext_inf, Ext_sup)
  pext[[5]]<-pext[[5]]+labs( x = expression(bar(epsilon))) + fundo_branco_grid + ylim(0, 1) +  ggtitle("d)")
  pext[[6]]<-pext[[6]]+labs( x = expression(bar(epsilon))) + fundo_branco_grid +  ggtitle("e)") 
  pext[[7]]<-pext[[7]]+labs( x = expression(bar(epsilon))) + fundo_branco_grid + ylim(0, 2000) +  ggtitle("f)")
  pcomp<-plot_function(Mean_comp, Inf_comp, Sup_comp)
  pcomp[[5]]<-pcomp[[5]]+labs( x = " ") + fundo_branco_grid + ylim(0, 1) +  ggtitle("a)")
  pcomp[[6]]<-pcomp[[6]]+labs( x = " ") + fundo_branco_grid +   ggtitle("b)")
  pcomp[[7]]<-pcomp[[7]]+labs( x = " ") + fundo_branco_grid + ylim(0, 2000) +  ggtitle("c)")
  plots <- list(pcomp[[5]],pcomp[[6]],pcomp[[7]],
                pext[[5]],pext[[6]],pext[[7]])
  
  legend<-pext[[1]]+theme(legend.position = "bottom",
                          legend.title = element_text(size=20),
                          legend.text  = element_text(size=20))
  shared_legend<-get_legend(legend)
  
  
  grobs <- list()
  widths <- list()
  
  for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  
  maxwidth <- do.call(grid::unit.pmax, widths)
  
  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  dev.off()
  
  #setwd("~/Documents")
  #pdf("fig4.pdf",width = 10,height=6.5)
  p <-grid.arrange(arrangeGrob(
    grobs[[1]],grobs[[2]],grobs[[3]],
    grobs[[4]],grobs[[5]],grobs[[6]],nrow = 2, ncol = 3),
    shared_legend, nrow=2, heights = c(10, 1))
  #dev.off()
  
}
fig5<-function(plot1, plot2){
  
  
  m<-0.002; p<-200;n<-1
  l1<-80
  objT1=read.table(paste0("Pop",p,"_Isol",l1,"_Mig",m,"_Temporal",n,"_.txt"),header=TRUE)
  l2<-20
  objT2=read.table(paste0("Pop",p,"_Isol",l2,"_Mig",m,"_Temporal",n,"_.txt"),header=TRUE)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  cores = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  ncor=length(cores)
  
  par(fig=c(0.025,0.5,0.40,1),new=T)
  
  
  
  p1<-plot(plot1[,1],plot1[,2], type="l",lwd=2, cex.main=1.5,
           main=expression("t"[h]%~~%"80%"),cex.lab=1.5,xaxt="n",xlab="",ylab="Balance (J)",ylim=c(0.0,1.0))
  barrier(objT1,0.25)
  lines(plot1[,1],plot1[,3],lty=3)
  lines(plot1[,1],plot1[,4],lty=3)
  p1<-rect(1980, -5, 2000, 5, density = NULL, angle = 45,
           col= rgb(0.502,0.502,0.502,alpha=0.25),border = NA, lty = par("lty"), lwd = par("lwd"))
  
  
  par(fig=c(0.48,0.9575,0.4,1),new=T)
  p3<-plot(plot2[,1],plot2[,2], type="l",lwd=2,cex.main=1.5,
           main=expression("t"[h]%~~%"20%"),cex.lab=1.5,xaxt="n",xlab="",ylab="",ylim=c(0.0,1.0))
  barrier(objT2,0.25)
  lines(plot2[,1],plot2[,3],lty=3)
  lines(plot2[,1],plot2[,4],lty=3)
  p3<-rect(1990, -5, 2000, 5, density = NULL, angle = 45,
           col= rgb(0.502,0.502,0.502,alpha=0.25),border = NA, lty = par("lty"), lwd = par("lwd"))
  
  
  
  par(fig=c(0.025,0.5,0.0,0.6),new=T)
  p2<-plot(plot1[,1],plot1[,5], type="l", lwd=2, cex.lab=1.5,xlab="Iterations", ylab=expression(alpha~" value"), ylim=c(-0.8,1.8))
  barrier(objT1,0.25)
  lines(plot1[,1],plot1[,6],lty=3)
  lines(plot1[,1],plot1[,7],lty=3)
  abline(h=c(0,1), col=c("red","red"),lty=c(2,2))
  p2<-rect(1980, -5, 2000, 5, density = NULL, angle = 45,
           col= rgb(0.502,0.502,0.502,alpha=0.25),border = NA, lty = par("lty"), lwd = par("lwd"))
  
  
  
  par(fig=c(0.48,0.9575,0.0,0.6),new=T)
  p4<-plot(plot2[,1],plot2[,5], type="l", lwd=2,cex.lab=1.5, xlab="Iterations", ylab="",ylim=c(-4,4.0))
  barrier(objT2,0.25)
  lines(plot2[,1],plot2[,6],lty=3)
  lines(plot2[,1],plot2[,7],lty=3)
  abline(h=c(0,1), col=c("red","red"),lty=c(2,2))
  p4<-rect(1990, -5, 2000, 5, density = NULL, angle = 45,
           col= rgb(0.502,0.502,0.502,alpha=0.25),border = NA, lty = par("lty"), lwd = par("lwd"))
  
  
  
}
fig6<-function(correlation_values){
  
  col2 <- colorRampPalette(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                      "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", 
                                      "#2166AC"))(200)  
                                      
  colnames(correlation_values)=c("mig","level","J","Age","Alpha","Asymr","Beta","Richness","Speciation")
  par(mfrow=c(2,2))
  layout(matrix(c(1,2,3,4,5,3), 2, 3, byrow = TRUE),widths=c(3,3,1))
  values<-correlation_values
  m<-c("0.005","0.01","0.0175","0.02","0.025","0.03","0.035","0.04","0.045","0.05","0.0575","0.06","0.065","0.07","0.08")
  dat<-subset(values,level %in% 0)
  dat<-dat[, -1]
  dat<-dat[, -1]
  dat<-dat[!rowSums(is.na(dat)),drop=F, ]
  corr <- round(cor(dat), 1)
  p.mat <- cor_pmat(corr)
  p0<-corrplot(corr, type = "full"
               ,method="color"
               ,addCoef.col = "black"
                 ,col=rev(col2)
               ,tl.col="black"
                 ,tl.srt=45
               ,cl.pos="n"
               ,tl.cex=1.2
               ,number.cex=1.5
               #,p.mat=p.mat
  ) 
  mtext(expression("t"[h]~"= 0"), at=1.25, line=0.0, cex=1.6)
  
  values<-correlation_values
  dat<-subset(values,level %in% 20)
  dat<-dat[, -1]
  dat<-dat[, -1]
  dat<-dat[!rowSums(is.na(dat)),drop=F, ]
  corr <- round(cor(dat), 1)
  p.mat <- cor_pmat(corr)
  p20<-corrplot(corr, type = "full"
                ,method="color"
                ,addCoef.col = "black"
                  ,col=rev(col2)
                ,tl.col="black"
                  ,tl.srt=45
                ,cl.pos="n"
                ,tl.cex=1.2
                ,number.cex=1.5
                #,p.mat=p.mat
  ) 
  mtext(expression("t"[h]%~~%"20%"), at=1.75, line=0.0, cex=1.6)
  
  plot.new()
  colorlegend(rev(col2), -5:5/5, align = 'l', cex = 1.8, xlim = c(0.0, 0.6),
              ylim = c(0.05, 0.95), vertical = T,main="corr")
  
  values<-correlation_values
  dat<-subset(values,level %in% 40)
  dat<-dat[, -1]
  dat<-dat[, -1]
  dat<-dat[!rowSums(is.na(dat)),drop=F, ]
  corr <- round(cor(dat), 1)
  p.mat <- cor_pmat(corr)
  p40<-corrplot(corr, type = "full"
                ,method="color"
                ,addCoef.col = "black"
                  ,col=rev(col2)
                ,tl.col="black"
                  ,tl.srt=45
                ,cl.pos="n"
                ,tl.cex=1.2
                ,number.cex=1.5
                #,p.mat=p.mat
  ) 
  mtext(expression("t"[h]%~~%"40%"), at=1.75, line=0.0, cex=1.6)
  
  values<-correlation_values
  dat<-subset(values,level %in% 80)
  dat<-dat[, -1]
  dat<-dat[, -1]
  dat<-dat[!rowSums(is.na(dat)),drop=F, ]
  corr <- round(cor(dat), 1)
  p.mat <- cor_pmat(corr)
  p80<-corrplot(corr, type = "full"
                ,method="color"
                ,addCoef.col = "black"
                  ,col=rev(col2)
                ,tl.col="black"
                  ,tl.srt=45
                ,cl.pos="n"
                ,tl.cex=1.2
                ,number.cex=1.5
                #,p.mat=p.mat
  )
  mtext(expression("t"[h]%~~%"80%"), at=1.75, line=0.0, cex=1.6)
  
}
barrier= function(objT,transparency){
  
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  cores = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  ncor=length(cores)
  
  tmax<-2000; tcorte<-0
  objm=objT[objT$tem<=tmax,]
  objm$nesp=objm$esp
  objm$abund1=(objm$abund1/objm$pop)*0.8/0.5
  objm$abund2=(objm$abund2/objm$pop)*0.8/0.5
  tempos=unique(objm$tem)
  espt=max(unique(objm$nesp))
  times=ceiling(espt/ncor)
  cor=rep(cores, times)
  migsv=length(unique(objm$mig))
  mig_min=min(unique(objm$mig))
  
  r<-0
  s<-10
  for(t in tempos){
    
    mig=unique(objm$mig[objm$tem==t])
    
    if(migsv>1){
      
      if(mig==mig_min){
        
        if(r==0){ r<-r+1; ri<-t;}else{rf<-t;s<-0}
        
      }else{
        r<-0
        
        if(s==0){
          #dev.off()
          p1<-rect(rf, -5, ri, 5, density = NULL, angle = 45,
                   col= rgb(0.502,0.502,0.502,alpha=transparency),border = NA, lty = par("lty"), lwd = par("lwd"))
        }
        s<-1
      }
    }
    
  }
}
plot_function<-function(Mean_list, Inf_list, Sup_list){
  
  i<-0
  pp<-list() 
  
  #myList <- list(total_spec,rich,Beta,Asym,J1,Alpha,Age)
  yy<-c("Speciation events",expression("Richness (N"[T]~")"),
        expression(beta~"diversity"),expression("Asymmetry ("*Delta*N*")"),
        "Balance (J)", expression(alpha~"value"),"Phylogeny Age") 
  
  
  i<-i+1
  #Mean_list <- list(sac,alp,div,spe,ext,J,Beta,MDCA,TajD,Het,D0,D2) 
  df1 <- melt(as.matrix(Mean_list[[1]]))
  colnames(df1) <- c("level", "mig", "value") 
  dfI1<- melt(as.matrix(Inf_list[[1]]))
  colnames(dfI1) <- c("level", "mig", "value") 
  dfS1 <- melt(as.matrix(Sup_list[[1]]))
  colnames(dfS1) <- c("level", "mig", "value") 
  
  #E_mig<-mig*ts
  pp[[i]]<-ggplot(df1, aes(x=mig,y=value)) +
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),# hjust = 0.5),
          axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+#geom_vline(xintercept=-level,alpha=0.6,col="gray")+
    geom_line(aes(colour=factor(level)),size=1)+
    geom_ribbon(aes(ymin=dfI1$value, ymax=dfS1$value,fill=factor(level)),alpha=0.1)+
    #geom_errorbar(aes(ymin=value-dfD$value, ymax=value+dfD$value,color=factor(level)), width=.9)+
    labs(fill = "level", colour= "level", x ="E[m]", y = yy[i])+
    #labs(title = "Same expected number of migrants")+#, subtitle = "B=2000 g=0.9", caption = "[E/m*t]")+
    scale_color_viridis_d(name=expression("Isolation time t"[h]), labels = c("0", expression(""%~~%"20%"), expression(""%~~%"40%"), expression(""%~~%"80%")),aesthetics = c("color", "fill"), option = 'viridis')
  
  
  
  
  
  i<-i+1
  #Mean_list <- list(sac,alp,div,spe,ext,J,Beta,MDCA,TajD,Het,D0,D2) 
  df2 <- melt(as.matrix(Mean_list[[2]]))
  colnames(df2) <- c("level", "mig", "value") 
  dfI2<- melt(as.matrix(Inf_list[[2]]))
  colnames(dfI2) <- c("level", "mig", "value") 
  dfS2 <- melt(as.matrix(Sup_list[[2]]))
  colnames(dfS2) <- c("level", "mig", "value") 
  
  pp[[i]]<-ggplot(df2, aes(x=mig,y=value)) +
    #geom_tile(aes(level, mig, fill = value))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),# hjust = 0.5),
          axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+#geom_vline(xintercept=-level,alpha=0.6,col="gray")+
    geom_line(aes(colour=factor(level)),size=1)+
    geom_ribbon(aes(ymin=dfI2$value, ymax=dfS2$value,fill=factor(level)),alpha=0.1)+
    #geom_errorbar(aes(ymin=value-dfD$value, ymax=value+dfD$value,color=factor(level)), width=.9)+
    labs(fill = "level", colour= "level", x ="E[m]", y = yy[i])+
    #labs(title = "Same expected number of migrants")+#, subtitle = "B=2000 g=0.9", caption = "[E/m*t]")+
    scale_color_viridis_d("Seabed depth (h)",aesthetics = c("color", "fill"), option = 'viridis')
  
  
  i<-i+1
  #Mean_list <- list(sac,alp,div,spe,ext,J,Beta,MDCA,TajD,Het,D0,D2) 
  df3 <- melt(as.matrix(Mean_list[[3]]))
  colnames(df3) <- c("level", "mig", "value") 
  dfI3<- melt(as.matrix(Inf_list[[3]]))
  colnames(dfI3) <- c("level", "mig", "value") 
  dfS3 <- melt(as.matrix(Sup_list[[3]]))
  colnames(dfS3) <- c("level", "mig", "value") 
  
  pp[[i]]<-ggplot(df3, aes(x=mig,y=value)) +
    #geom_tile(aes(level, mig, fill = value))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),# hjust = 0.5),
          axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+#geom_vline(xintercept=-level,alpha=0.6,col="gray")+
    geom_line(aes(colour=factor(level)),size=1)+
    geom_ribbon(aes(ymin=dfI3$value, ymax=dfS3$value,fill=factor(level)),alpha=0.1)+
    #geom_errorbar(aes(ymin=value-dfD$value, ymax=value+dfD$value,color=factor(level)), width=.9)+
    labs(fill = "level", colour= "level", x ="E[m]", y = yy[i])+
    #labs(title = "Same expected number of migrants")+#, subtitle = "B=2000 g=0.9", caption = "[E/m*t]")+
    scale_color_viridis_d("Seabed depth (h)",aesthetics = c("color", "fill"), option = 'viridis')
  
  i<-i+1
  #Mean_list <- list(sac,alp,div,spe,ext,J,Beta,MDCA,TajD,Het,D0,D2) 
  df4 <- melt(as.matrix(Mean_list[[4]]))
  colnames(df4) <- c("level", "mig", "value") 
  dfI4<- melt(as.matrix(Inf_list[[4]]))
  colnames(dfI4) <- c("level", "mig", "value") 
  dfS4 <- melt(as.matrix(Sup_list[[4]]))
  colnames(dfS4) <- c("level", "mig", "value") 
  
  pp[[i]]<-ggplot(df4, aes(x=mig,y=value)) +
    #geom_tile(aes(level, mig, fill = value))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),# hjust = 0.5),
          axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+#geom_vline(xintercept=-level,alpha=0.6,col="gray")+
    geom_line(aes(colour=factor(level)),size=1)+
    geom_ribbon(aes(ymin=dfI4$value, ymax=dfS4$value,fill=factor(level)),alpha=0.1)+
    #geom_errorbar(aes(ymin=value-dfD$value, ymax=value+dfD$value,color=factor(level)), width=.9)+
    labs(fill = "level", colour= "level", x ="E[m]", y = yy[i])+
    #labs(title = "Same expected number of migrants")+#, subtitle = "B=2000 g=0.9", caption = "[E/m*t]")+
    scale_color_viridis_d("Seabed depth (h)",aesthetics = c("color", "fill"), option = 'viridis')
  
  i<-i+1
  #Mean_list <- list(sac,alp,div,spe,ext,J,Beta,MDCA,TajD,Het,D0,D2) 
  df5 <- melt(as.matrix(Mean_list[[5]]))
  colnames(df5) <- c("level", "mig", "value") 
  dfI5<- melt(as.matrix(Inf_list[[5]]))
  colnames(dfI5) <- c("level", "mig", "value") 
  dfS5 <- melt(as.matrix(Sup_list[[5]]))
  colnames(dfS5) <- c("level", "mig", "value") 
  
  pp[[i]]<-ggplot(df5, aes(x=mig,y=value)) +
    #geom_tile(aes(level, mig, fill = value))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),# hjust = 0.5),
          axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+#geom_vline(xintercept=-level,alpha=0.6,col="gray")+
    geom_line(aes(colour=factor(level)),size=1)+
    geom_ribbon(aes(ymin=dfI5$value, ymax=dfS5$value,fill=factor(level)),alpha=0.1)+
    #geom_errorbar(aes(ymin=value-dfD$value, ymax=value+dfD$value,color=factor(level)), width=.9)+
    labs(fill = "level", colour= "level", x ="E[m]", y = yy[i])+
    #labs(title = "Same expected number of migrants")+#, subtitle = "B=2000 g=0.9", caption = "[E/m*t]")+
    scale_color_viridis_d("Seabed depth (h)",aesthetics = c("color", "fill"), option = 'viridis')
  
  i<-i+1
  #Mean_list <- list(sac,alp,div,spe,ext,J,Beta,MDCA,TajD,Het,D0,D2) 
  df6 <- melt(as.matrix(Mean_list[[6]]))
  colnames(df6) <- c("level", "mig", "value") 
  dfI6<- melt(as.matrix(Inf_list[[6]]))
  colnames(dfI6) <- c("level", "mig", "value") 
  dfS6 <- melt(as.matrix(Sup_list[[6]]))
  colnames(dfS6) <- c("level", "mig", "value") 
  
  pp[[i]]<-ggplot(df6, aes(x=mig,y=value)) +
    #geom_tile(aes(level, mig, fill = value))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),# hjust = 0.5),
          axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+#geom_vline(xintercept=-level,alpha=0.6,col="gray")+
    geom_line(aes(colour=factor(level)),size=1)+ geom_hline(yintercept=1, linetype="dashed", color = "red")+
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    geom_ribbon(aes(ymin=dfI6$value, ymax=dfS6$value,fill=factor(level)),alpha=0.1)+
    #geom_errorbar(aes(ymin=value-dfD$value, ymax=value+dfD$value,color=factor(level)), width=.9)+
    labs(fill = "level", colour= "level", x ="E[m]", y = yy[i])+
    #labs(title = "Same expected number of migrants")+#, subtitle = "B=2000 g=0.9", caption = "[E/m*t]")+
    scale_color_viridis_d("Seabed depth (h)",aesthetics = c("color", "fill"), option = 'viridis')
  
  i<-i+1
  #Mean_list <- list(sac,alp,div,spe,ext,J,Beta,MDCA,TajD,Het,D0,D2) 
  df7 <- melt(as.matrix(Mean_list[[7]]))
  colnames(df7) <- c("level", "mig", "value") 
  dfI7<- melt(as.matrix(Inf_list[[7]]))
  colnames(dfI7) <- c("level", "mig", "value") 
  dfS7 <- melt(as.matrix(Sup_list[[7]]))
  colnames(dfS7) <- c("level", "mig", "value") 
  
  pp[[i]]<-ggplot(df7, aes(x=mig,y=value)) +
    #geom_tile(aes(level, mig, fill = value))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),# hjust = 0.5),
          axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+#geom_vline(xintercept=-level,alpha=0.6,col="gray")+
    geom_line(aes(colour=factor(level)),size=1)+
    geom_ribbon(aes(ymin=dfI7$value, ymax=dfS7$value,fill=factor(level)),alpha=0.1)+
    #geom_errorbar(aes(ymin=value-dfD$value, ymax=value+dfD$value,color=factor(level)), width=.9)+
    labs(fill = "level", colour= "level", x ="E[m]", y = yy[i])+
    #labs(title = "Same expected number of migrants")+#, subtitle = "B=2000 g=0.9", caption = "[E/m*t]")+
    scale_color_viridis_d("Seabed depth (h)",aesthetics = c("color", "fill"), option = 'viridis')
  
  
  return(pp)
}
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
fundo_branco_grid=theme(panel.grid.major = element_line(colour="gray")
                        , panel.grid.minor =element_blank()
                        , panel.background = element_blank()
                        , axis.line = element_line(colour = "black")
                        , panel.border = element_rect(colour = "black", fill=NA, size=1)
                        , plot.title = element_text(hjust = 0.005, size=15)
                        , plot.tag = element_text(size=20)
                        , plot.tag.position = c(0.935, 0.90)
                        , axis.title.x = element_text(size = 30)
                        , axis.title.y = element_text(size = 20)
                        , legend.text = element_text(size = 20)
)
################################################################################


# UTILIZING THE PROVIDED DATA TO GENERATE THE FIGURES #### 

#Set the directory to load the .RData files below
setwd("~/Data/Results")
# COMPLETE TREE DATA
Comp_mean<-readRDS("Comp_mean.RData")
Comp_inf<-readRDS("Comp_inf.RData")
Comp_sup<-readRDS("Comp_sup.RData")

# EXTANT TREE DATA
Ext_mean<-readRDS("Ext_mean.RData")
Ext_inf<-readRDS("Ext_inf.RData")
Ext_sup<-readRDS("Ext_sup.RData")

# J and Alpha over evolutionary time
plot80<-readRDS("plot80.RData")
plot20<-readRDS("plot20.RData")

# CORRELATION
correlation_values<-readRDS("correlation_values.RData")



#### Figures from Freitas et al. (2023) #####
fig3()
dev.off()

fig4()
dev.off()

fig6(correlation_values)
dev.off()

fig5(plot80, plot20)
dev.off()


##############################################






# Heavly Modified in: 04/05/23 ####
