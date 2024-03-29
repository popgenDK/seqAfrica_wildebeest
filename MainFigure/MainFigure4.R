## only for Loter output

require(data.table)
require(ggplot2)
require(tidyverse)

setwd('/Users/pls394/Documents/Project/Wildebeest/Loter/plot_sanity_check/plot_chr_geno_correction_nyassa')
loter_list = read.table("/Users/pls394/Documents/Project/Wildebeest/Loter/plot_sanity_check/nyassa_list.out",header=F) # nyassa
loter_list = read.table("/Users/pls394/Documents/Project/Wildebeest/Loter/plot_sanity_check/loter_list.out",header=F) # brindled
pos_list = read.table("/Users/pls394/Documents/Project/Wildebeest/Loter/plot_sanity_check/pos_list.out",header=F)

pos_list2 = as.character(pos_list$V1)
loter_list2 = as.character(loter_list$V1)
t = do.call(rbind, lapply(pos_list2, function(x) read.table(x, header=F)))
s= do.call(rbind,lapply(loter_list2,function(x) t(fread(x,sep=" ",header=F)) ))
cat("Reading input done!\n")
dat = cbind(t,s)
colnames(dat)= c('chr','pos',paste0('hap',1:10))


number = apply(s,1,function(x) {length(which(x==1))})
unique(t[which(number>=10),1])

### Now go through every ind ###

summary <- data.frame(ind=numeric(),hap=numeric(),ancestry=numeric(),distance=numeric())
summary.gt <- data.frame(ind=numeric(),ancestry=numeric(),distance=numeric())

for (num in 1:5 ) { # go  through ind
    
    ind.index = num*2 + 1 ## ind index
    cat("Individual-",num,"\n")
    
    y.start = 0
    y.space = 1
    rect.length = 0.4 #rectangle length/2
    
    triangle.d = data.frame(x1 = numeric(), # for homo of blue ancestry
                            x2 = numeric(),
                            y1 = numeric(),
                            y2 = numeric()
    )
    
    segment.d = data.frame(x=numeric(), # for hetero of black ancestry
                           y=numeric(),
                           xend=numeric(),
                           yend=numeric()
    )
    
    segment.d.2 = data.frame(x=numeric(), # for heomo of black ancestry
                           y=numeric(),
                           xend=numeric(),
                           yend=numeric()
    )
    
    
    for (chr in unique(dat$chr)) { # go through chr
        chr.index = which(dat$chr == chr)
        
        
        hap1  = as.numeric( dat[chr.index,ind.index])
        hap2 = as.numeric( dat[chr.index,ind.index+1])
        
        ### deal with hap1 ###
        rle_hap1 = rle(hap1)
        end_hap1 = cumsum(rle_hap1$lengths)
        start_hap1 = c(1, lag(end_hap1)[-1] + 1)
        type_hap1 = rle_hap1$values
        summary.hap1 <- data.frame(ind=num,hap=1,ancestry=type_hap1,
                                   distance=dat[chr.index,'pos'][end_hap1]-dat[chr.index,'pos'][start_hap1])
        summary <- rbind(summary,summary.hap1)
        ### end of hap1 ###
        
        ### deal with hap2 ###
        rle_hap2 = rle(hap2)
        end_hap2 = cumsum(rle_hap2$lengths)
        start_hap2 = c(1, lag(end_hap2)[-1] + 1)
        type_hap2 = rle_hap2$values
        summary.hap2 <- data.frame(ind=num,hap=2,ancestry=type_hap2,
                                   distance=dat[chr.index,'pos'][end_hap2]-dat[chr.index,'pos'][start_hap2])
        summary <- rbind(summary,summary.hap2)
        ### end of hap2 ###
        
        ### this part deal with genotype, hap1+ hap2 ####
        geno = hap1 + hap2
        #chunk.geno = split(geno, cumsum(c(1, diff(geno) != 0)))
        # df1, data.frame(pos,gt=0,1,1,2,0)
        df1 = data.frame(position= dat[chr.index,"pos"],gt=geno)
        rle_x = rle(geno)
        end = cumsum(rle_x$lengths)
        start = c(1, lag(end)[-1] + 1)
        #type =as.numeric(sapply(chunk.geno,"[[",1))
        type = rle_x$values
        
        # df2, data.frame(start,end, type=0,1,2)
        df2 = data.frame(start=dat[chr.index,'pos'][start], end=dat[chr.index,'pos'][end],type=type)
        ##### end of genotype ####
        
        
        #### this estmiate chunks 0,1 ####
        
        # here, we are mainly interested in minor ancestry (black), so either gt==1 (hetro) or gt ==2 (homo of black)
        df1.binary = df1
        df1.binary[which(df1$gt ==2),'gt'] =1
        
        rle_x.binary = rle(df1.binary$gt)
        end.binary = cumsum(rle_x.binary$lengths)
        start.binary = c(1, lag(end.binary)[-1] + 1)
        #type =as.numeric(sapply(chunk.geno,"[[",1))
        type.binary = rle_x.binary$values
        
        df2.binary = data.frame(start=dat[chr.index,'pos'][start.binary], end=dat[chr.index,'pos'][end.binary],
                                type=type.binary)
        df2.binary$lenght = df2.binary$end - df2.binary$start
        assign(paste(chr,"df2.binary",sep="."),df2.binary)
        summary.gt.binary <- data.frame(ind=num,ancestry=df2.binary$type,
                                   distance=df2.binary$end - df2.binary$start)
        summary.gt <- rbind(summary.gt,summary.gt.binary)
        
       
        
        ### for plot ###
        triangle.new = data.frame(x1 = df2[df2$type==0,'start'],
                                  x2 = df2[df2$type==0,"end"],
                                  y1 = y.start-rect.length,
                                  y2 = y.start +rect.length
        )
        triangle.d = rbind(triangle.d, triangle.new)
        
        if(length(which(df1$gt=="1")) >=1 ) { # hetero for black and blue
          segment.new = data.frame(x=df1[df1$gt=="1","position"],
                                 y=y.start-rect.length,
                                 xend =df1[df1$gt=="1","position"],
                                 yend=y.start+rect.length
          )
          segment.d = rbind(segment.d,segment.new)
        } 
        
        if(length(which(df1$gt=="2")) >=1 ) { # homo for black
          segment.new.2 = data.frame(x=df1[df1$gt=="2","position"],
                                   y=y.start-rect.length,
                                   xend =df1[df1$gt=="2","position"],
                                   yend=y.start+rect.length
          )
          segment.d.2 = rbind(segment.d.2,segment.new.2)
        }
        
        y.start = y.start + y.space
    } # end of loop of chr   

    #### plot here #####
    # g1 = ggplot() + 
    #     scale_x_continuous(name="x") + 
    #     scale_y_continuous(name="y") +
    #     geom_rect(data=triangle.d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="lightgray", fill=alpha("lightgray", 1))
    # 
    # g2 = g1 + geom_segment(data=segment.d.2, mapping= aes(x=x, y =y, xend =xend,yend=yend,col="black"),size = 0.015)
    # 
    # g3 = g2 + geom_segment(data=segment.d, mapping= aes(x=x, y =y, xend =xend,yend=yend),col="magenta2",size = 0.015)
    # 

    # 
    #     
    # bitmap(paste0("ind-",num,".GT.ancestry.correciton.new.png"),h=4,w=5,res=400)
    # g4 =g3 + myblanktheme +  scale_y_continuous(name = "Chromosomes", # remove y title
    #                                        breaks = 0:(y.start-1)) +
    #                     scale_x_continuous(name=NULL,
    #                                        breaks = seq(0,180000000,20000000),
    #                                        labels = paste0( seq(0,180,20), c(rep("",9),"Mb"))) 
    
    
    #g5 = g4 +   scale_color_manual(name='Genotype',
    #                               breaks=c('Homozygote of blue wildebeest allele', 'Homogzygote of black wildbeest allele', 'Heterozygote'),
    #                               values=c('Homozygote of blue wildebeest allele'='lightgray', 'Homogzygote of black wildbeest allele'='black', 'Heterozygote'='magenta2'))                    
    #print(g5)
    #dev.off()
    
    scale_custom <- list(
        scale_color_manual(breaks=c( 'Homogzygote of black wildbeest allele', 'Heterozygote'),
                           values=c( 'Homogzygote of black wildbeest allele'='black', 'Heterozygote'='magenta2')),
        scale_fill_manual(breaks=c('Homozygote of blue wildebeest allele'),values=c('Homozygote of blue wildebeest allele'='lightgray'))
    )
    
    
     # ggplot() + 
     #    scale_x_continuous(name="x") + 
     #    scale_y_continuous(name="y") +
     #    geom_rect(data=triangle.d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill="Homozygote of blue wildebeest allele"),color="lightgray") +
     #    scale_fill_manual(name='Genotype',breaks=c('Homozygote of blue wildebeest allele'),values=c('Homozygote of blue wildebeest allele'='lightgray')) +
     #    geom_segment(data=segment.d.2, mapping= aes(x=x, y =y, xend =xend,yend=yend,color="Homogzygote of black wildbeest allele"),size = 0.015) +
     #    geom_segment(data=segment.d, mapping= aes(x=x, y =y, xend =xend,yend=yend,color="Heterozygote"),size = 0.015)  +
     #    scale_color_manual(name='Genotype',
     #                        breaks=c( 'Homogzygote of black wildbeest allele', 'Heterozygote'),
     #                        values=c( 'Homogzygote of black wildbeest allele'='black', 'Heterozygote'='magenta2'))  +
     #    myblanktheme +  scale_y_continuous(name = "Chromosomes", # remove y title
     #                                        breaks = 0:(y.start-1)) +
     #    scale_x_continuous(name=NULL,
     #                        breaks = seq(0,180000000,20000000),
     #                        labels = paste0( seq(0,180,20), c(rep("",9),"Mb"))) + theme(legend.position=c(.8,.5),legend.background=element_blank())
     
     
     # ggplot() + 
     #     scale_x_continuous(name="x") + 
     #     scale_y_continuous(name="y") +
     #     geom_rect(data=triangle.d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill="Homozygote of blue wildebeest allele"),color="lightgray") +
     #   
     #     geom_segment(data=segment.d.2, mapping= aes(x=x, y =y, xend =xend,yend=yend,color="Homogzygote of black wildbeest allele"),size = 0.015) +
     #     geom_segment(data=segment.d, mapping= aes(x=x, y =y, xend =xend,yend=yend,color="Heterozygote"),size = 0.015)  +
     #     
     #     myblanktheme +  scale_y_continuous(name = "Chromosomes", # remove y title
     #                                        breaks = 0:(y.start-1)) +
     #     scale_custom +
     #     scale_x_continuous(name=NULL,
     #                        breaks = seq(0,180000000,20000000),
     #                        labels = paste0( seq(0,180,20), c(rep("",9),"Mb"))) + theme(legend.position=c(.8,.7)) + theme(legend.key=element_blank())
     
     bitmap(paste0("/Users/pls394/Documents/Project/Wildebeest/Loter/MS_figure/","ind-",num,".GT.ancestry.correciton.new.png"),h=4,w=6,res=400)
     
     
     myblanktheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "white"),
                          axis.title=element_text(size=14,face="bold"), axis.ticks.y=element_blank())
     
     fig = ggplot() + 
         scale_x_continuous(name="x") + 
         scale_y_continuous(name="y") +
         geom_rect(data=triangle.d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill="Homozygote of blue wildebeest allele"),color="lightgray") +
         
         geom_segment(data=segment.d.2, mapping= aes(x=x, y =y, xend =xend,yend=yend,color="Homogzygote of black wildbeest allele"),size = 0.01) +
         geom_segment(data=segment.d, mapping= aes(x=x, y =y, xend =xend,yend=yend,color="Heterozygote"),size = 0.01)  +
         
         myblanktheme +  scale_y_continuous(name = "Chromosomes", # remove y title
                                            breaks = 0:(y.start-1),labels=c(1,2,4:29)) +
         scale_custom +
         scale_x_continuous(name=NULL,
                            breaks = seq(0,180000000,20000000),
                            labels = paste0( seq(0,180,20), c(rep("",9),"Mb"))) + theme(legend.position=c(.8,.7)) + theme(legend.position = 'none')
    plot(fig)                                                                                                                         
    dev.off()
    
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("topleft", legend =c('homozygote of black wildebeest allele', 'homozygote of blue wildebeest allele', 'heterozygotes'), pch=15, pt.cex=3, cex=1.5, bty='n',
           col = c('black', 'lightgray', 'magenta2'))
    mtext("Genotypes", at=0.2, cex=2)
    ### end of plot
}
    


########## end of parsing data #########
 



### plot summary, genotyp
par(mfrow=c(5,4),mai = c(0.3, 0.3, 0.1, 0.3))
for (ind in 1:5) {
    for (anc in 0:1) {
        index = which(summary.gt$ind==ind &  summary.gt$ancestry==anc)
        x = as.numeric(summary.gt[index,3])/1000000
        hist(x,main="",xlab="",breaks=50)
            #lines( density(x),col="red" )
        boxplot(x)
        print(max(x,na.rm = T))
    }
}

k
