




getAnc <- function(x,df,rev=FALSE){

    anc <-  df[df$to==tail(x,1),"from"]
    if(length(anc)>2){
        cat("2> anc. To",tail(x,1),"\n")
        print(anc)
    }
    if(length(anc)==2){
        if(rev)
            return( c(getAnc(c(x,anc[1]),df=df),getAnc(anc[2],df=df)))
        else
            return( c(getAnc(c(x,anc[2]),df=df),getAnc(anc[1],df=df)))
    }
    if(length(anc)==0)
        return(x)
  
    getAnc(c(x,anc),df=df)  
}

# returns the node with time fuckup

timepolice <- function(df){

    nodes <- unique(df$to)
    ancestors <-lapply(nodes,getAnc,df=df)
    ancestors2 <-lapply(nodes,getAnc,df=df,rev=T)

    ##is the first direct ancestor node also part of the path of the second ancestror
    timeFuck <- sapply(ancestors,function(x) x[2]%in%x[-2] )
    ## same but swich ancestors
    timeFuck2 <- sapply(ancestors2,function(x) x[2]%in%x[-2] )

    nodes[timeFuck | timeFuck2]
}


if(FALSE){
 
    source("/home/albrecht/projects/admixer2022/scripts/timepolice.R")
    df <- read.table("/home/krishang/to_anders_timepolice.txt",head=T,as.is=T)
fed(df)
    
if(length( res<-timepolice(df) )>0 )
    cat("TIME VIOLATION with nodes:",res,"\n")


}
