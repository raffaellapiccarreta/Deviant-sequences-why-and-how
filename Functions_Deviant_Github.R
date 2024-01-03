

obtain.deviations<-function(Partition,Distance,NN=c(5,10,15,20),
                            ref.percentiles=c(0.5,0.75,0.8,0.85,0.9,0.95,0.99)){
  
  # VERIFYING CONSISTENCY OF INPUTS
  
  # stop if inputs not specified or not existing
  if(isTRUE(missing(Partition))){
    stop("Argument 'Partition' missing")} 
  if(isTRUE(missing(Distance))){
    stop("Argument 'Distance' missing")} 
  ge <- (ls(name=.GlobalEnv))
  name.part<-deparse1(substitute(Partition))
  # if(!(name.part %in% ge)){stop(paste0(name.part," does not exist"))} 
  name.dist<-deparse1(substitute(Distance))
  # if(!(name.dist %in% ge)){stop(paste0(name.dist," does not exist"))} 
  
  # check type of object
  if((is.matrix(Partition) | is.data.frame(Partition)) && ncol(Partition)>1){
    stop(paste0(name.part,"is a ",class(Partition),"with more than 1 column"))}
  if((is.matrix(Partition) | is.data.frame(Partition)) && ncol(Partition)==1){
    Partition<-as.factor(Partition[,1])}
  if(!is.atomic(Partition) & !is.factor(Partition)){
    stop(paste0("Argument Partition cannot be a ",class(Partition)))}
  if(is.atomic(Partition)){Partition<-as.factor(Partition)}
  
  if (inherits(Distance, "dist")){
    Distance<-as.matrix(Distance)
  } else if(is.matrix(Distance)){
    if(ncol(Distance)!=nrow(Distance)){
      stop("Argument Distance should be a squared matrix or a 'dist' object")
    }
  }
  if(sum(abs(Distance-t(Distance))) !=0){
    stop("Argument Distance should be a squared matrix or a 'dist' object")
  }
  if(nrow(Distance) != length(Partition)){
    stop("Argument Distance should be a squared matrix or a 'dist' object")
  }
  if(sum(is.na(Distance))>0  |  sum(is.na(Partition))>0){
    stop(name.part," or ",name.dist," have missing values")
  }
  if(length(levels(Partition))<2){
    stop(name.part," must have at least 2 levels")
  }
  
  min.in.cl<-min(table(Partition))-1
  if(sum(NN<=min.in.cl)==0){
    stop("At least one cluster has too small size for the specified NNs")
  }
  if(sum(NN<=min.in.cl)<length(NN)){
    cat("Warning: \nNNs cannot be higher than clusters' sizes \n Only admissible NNs are considered:\n",file=stderr())
    NN<-NN[NN<=min.in.cl]
    write.table(format(t(as.matrix(NN)), justify="left"),
                row.names=F, col.names=F, quote=F)  }
  
  # EVALUATION OF DISSIMILARITY USING DIFFERENT CRITERIA
  cl.labs<-levels(Partition)
  diag(Distance)<-NA
  
  # distance from own cluster
  # distance from cluster's medoid
  avedist.own<-avedist.others<-rep(NA,nrow(Distance))
  dist.medoid<-silhouette.mine<-rep(NA,nrow(Distance))
  for(k in 1:length(cl.labs)){
    sel.k<-Partition==cl.labs[k]
    dist.k<-Distance[sel.k,]
    avedist.own[sel.k]<-apply(Distance[sel.k,sel.k],1,mean,na.rm=T)
  }

  # distance from medoids
  dist.medoid<-rep(NA,length(Partition))
  diag(Distance)<-0
  for(k in 1:length(cl.labs)){
    usedist<-Distance[Partition==cl.labs[k],Partition==cl.labs[k]]
    dc <- disscenter(usedist,medoids.index="first")
    dist.medoid[Partition==cl.labs[k]]<-usedist[,dc]
  }
  
  # distance from Nearest Neighbours
  func.use<-function(x,nns){y=x[order(x)]; y[is.na(y)==F];  y<-max(y[1:nns])}
  
  distfromNN_inside<-NULL
  diag(Distance)<-NA
  for(nn in NN){
    inside.nn<-rep(NA,nrow(Distance))
    for(k in 1:length(cl.labs)){
      sel.k<-Partition==cl.labs[k]
      dist.k<-Distance[sel.k,sel.k]
      inside.nn[sel.k]<-sapply(as.data.frame(dist.k),function(x) func.use(x,nn))
    }
    distfromNN_inside<-cbind(distfromNN_inside,inside.nn)
  }
  colnames(distfromNN_inside)<-paste0("DistfromNN.",NN)
  
  Deviations<-data.frame(Cluster=Partition,Avedist.own=avedist.own,
                         Dist.FromMedoid=dist.medoid,distfromNN_inside)
  
  # CALCULATING REFERENCE MEASURES
  
  gen.stats<-apply(Deviations[,2:ncol(Deviations)],2,
                   function(x) c(mean(x),sd(x),quantile(x,probs=ref.percentiles)))
  gen.stats<-as.data.frame(gen.stats)
  gen.stats$Type<-"General"
  gen.stats$Summary<-c("Average","Std.Dev",paste0("Pct",ref.percentiles*100))
  gen.stats$Cluster<-NA
  rownames(gen.stats)<-paste0("Gen_",gen.stats$Summary)
  
  bycl.stats<-aggregate(Deviations[,2:ncol(Deviations)],
                        by=list(Cluster=Deviations$Cluster),"mean")
  bycl.stats$Type<-"ByCl"
  bycl.stats$Summary<-"Average"
  
  bycl.stats.s<-aggregate(Deviations[,2:ncol(Deviations)],
                          by=list(Cluster=Deviations$Cluster),"sd")
  bycl.stats.s$Type<-"ByCl"
  bycl.stats.s$Summary<-"Std.Dev"
  bycl.stats<-rbind(bycl.stats,bycl.stats.s)
  for(pct in ref.percentiles){
    bycl.stats.s<-aggregate(Deviations[,2:ncol(Deviations)],
                            by=list(Cluster=Deviations$Cluster),
                            "quantile",probs=pct)
    bycl.stats.s$Type<-"ByCl"
    bycl.stats.s$Summary<-paste0("Pct",pct*100)
    bycl.stats<-rbind(bycl.stats,bycl.stats.s)
  }
  rownames(bycl.stats)<-paste0("ByCl_",bycl.stats$Cluster,"_",bycl.stats$Summary)
  bycl.stats<-bycl.stats[,colnames(gen.stats)]
  
  All.Stats<-rbind(gen.stats,bycl.stats[,colnames(gen.stats)])
  
  
  # save results
  final.list<-list(Deviations,All.Stats,createdby="obtain.deviations")
  names(final.list)<-c("Deviations","Reference.Values","createdby")
  return(final.list)
}

obtain.rel_deviations<-function(info.deviations,criterion,
                        ref.location,ref.scale,ref.noise,
                        Distance){
  print.with.c<-function(msg,vec){
    for(k in 1:length(vec)){msg<-paste(msg,vec[k])}
    return(msg)
  }
  
  if(isTRUE(missing(info.deviations))){
    stop("Argument 'info.deviations' missing")} 
  if(!("createdby" %in% names(info.deviations))){
    stop("Argument 'info.deviations' should be the output of function obtain.deviations")} 
  if(isTRUE(missing(criterion))){
    stop("Argument 'criterion' missing")} 
  if((isFALSE(missing(ref.location)) | isFALSE(missing(ref.scale)) ) &
     (isFALSE(missing(ref.noise)) )){
    stop("Both (*location or *scale) and *noise references are specified")
  } 
  if(isTRUE(missing(ref.location)) & isTRUE(missing(ref.scale))  &
     isTRUE(missing(ref.noise)) ){
    stop("No references specified")
  } 
  
  if(isFALSE(missing(ref.location)) | isFALSE(missing(ref.scale))  ){
    vec.crit=colnames(info.deviations$Deviations)[-1]}
  if(isFALSE(missing(ref.noise))){
    vec.crit<-colnames(info.deviations$Deviations)
    vec.crit<-vec.crit[substring(vec.crit,1,10)=="DistfromNN"]}
  
  if(length(criterion)>1){
    stop(print.with.c(msg="Only one 'criterion' can be specified; options are\n", 
                      vec=vec.crit))  }
  if(!(criterion %in% vec.crit)){
    stop(print.with.c(msg="Wrong specification for 'criterion'; admissible options are\n", 
                      vec=vec.crit)) }
  
  ref.stats<-info.deviations$Reference.Values
  cases.values<-info.deviations$Deviations[,criterion]
  
  # center values if correctly required
  if(isFALSE(missing(ref.location))){
    if(length(ref.location)>1 | !is.character(ref.location)){
      stop("Argument 'ref.location' should be a 1-element character vector")} 
    
    ok.stats<-ref.stats$Summary[duplicated(ref.stats$Summary)==F]
    ok.stats<-ok.stats[ok.stats!="Std.Dev"]
    if(!(ref.location %in% c(ok.stats))){
      stop(print.with.c(msg="Wrong 'ref.location' specification; admissible options are:\n",
                        vec=ok.stats))
    }
    if(is.null(names(ref.location))){
      stop("Argument 'ref.location' should be a named 1-element character vector")} 
    if(!(names(ref.location) %in% c("General","ByCl"))){
      stop("Wrong name for 'ref.location': admissible options are:\n General ByCl")
    } 
    if(names(ref.location)=="General"){
      sel.loc<-ref.stats$Type=="General" & ref.stats$Summary==ref.location
      cases.values<-cases.values-ref.stats[sel.loc,criterion]
    }
    if(names(ref.location)=="ByCl"){
      sel.loc<-ref.stats$Type=="ByCl" & ref.stats$Summary==ref.location
      sel.loc<-ref.stats[sel.loc,]
      for(k in levels(info.deviations$Deviations$Cluster)){
        sel.k<-info.deviations$Deviations$Cluster==k
        cases.values[sel.k]<-cases.values[sel.k]-sel.loc[sel.loc$Cluster==k,criterion]
      }
    }
  } # closes ref.location
  
  # scale values if correctly required
  if(isFALSE(missing(ref.scale))){
    if(length(ref.scale)>1 | !is.character(ref.scale)){
      stop("Argument 'ref.scale' should be a 1-element character vector")} 
    
    ok.stats<-ref.stats$Summary[duplicated(ref.stats$Summary)==F]
    if(!(ref.scale %in% c(ok.stats))){
      stop(print.with.c(msg="Wrong 'ref.scale' specification; admissible options are:\n",
                        vec=ok.stats))
    }
    if(is.null(names(ref.scale))){
      stop("Argument 'ref.scale' should be a named 1-element character vector")} 
    if(!(names(ref.scale) %in% c("General","ByCl"))){
      stop("Wrong name for 'ref.scale': admissible options are:\n General ByCl")
    } 
    if(names(ref.scale)=="General"){
      sel.scale<-ref.stats$Type=="General" & ref.stats$Summary==ref.scale
      cases.values<-cases.values/ref.stats[sel.scale,criterion]
    }
    if(names(ref.scale)=="ByCl"){
      sel.scale<-ref.stats$Type=="ByCl" & ref.stats$Summary==ref.scale
      sel.scale<-ref.stats[sel.scale,]
      for(k in levels(info.deviations$Deviations$Cluster)){
        sel.k<-info.deviations$Deviations$Cluster==k
        cases.values[sel.k]<-cases.values[sel.k]/sel.scale[sel.scale$Cluster==k,criterion]
      }
    }
  } # closes ref.scale
  
  # find noisy cases if correctly required
  if(isFALSE(missing(ref.noise))){
    if(length(ref.noise)>1 | !is.character(ref.noise)){
      stop("Argument 'ref.noise' should be a 1-element character vector")} 
    
    ok.stats<-ref.stats$Summary[duplicated(ref.stats$Summary)==F]
    ok.stats<-ok.stats[!(ok.stats %in% c("Average","Std.Dev"))]
    if(!(ref.noise %in% c(ok.stats))){
      stop(print.with.c(msg="Wrong 'ref.noise' specification; admissible options are:\n",
                        vec=ok.stats))
    }
    if(isTRUE(missing(Distance))){
      stop("Argument 'Distance' is needed when 'ref.noise' is specified")} 
    if (inherits(Distance, "dist")){
      Distance<-as.matrix(Distance)
    } else if(is.matrix(Distance)){
      if(ncol(Distance)!=nrow(Distance)){
        stop("Argument Distance should be a squared matrix or a 'dist' object")
      }
    }
    if(sum(abs(Distance-t(Distance))) !=0){
      stop("Argument Distance should be a squared matrix or a 'dist' object")
    }
    if(nrow(Distance) != nrow(info.deviations$Deviations)){
      stop("'Distance' has a dimension incompatible with data in 'info.deviations'")
    }
    if(sum(is.na(Distance))>0){
      stop("'Distance'  has missing values")
    }
    
    if(!is.null(names(ref.noise)) && names(ref.noise)!="General"){
      cat("Warning: With 'ref.noise' the reference percentile is based on the entire sample",
          file=stderr())} 
    
    
    sel.noise<-ref.stats$Type=="General" & ref.stats$Summary==ref.noise
    sel.noise<-ref.stats[sel.noise,criterion]
    NN<-as.numeric(substring(criterion,12))
    TypeCase<-rep("Noise",length(cases.values))
    
    diag(Distance)<-NA
    # identify core points
    for(k in levels(info.deviations$Deviations$Cluster)){
      sel.k<-info.deviations$Deviations$Cluster==k
      dist.k<-Distance[sel.k,sel.k]<=sel.noise
      forcore_k<-apply(dist.k,1,sum,na.rm=T)
      TypeCase[sel.k][forcore_k>=NN]<-"Core"
    }
    # identify border points
    for(k in levels(info.deviations$Deviations$Cluster)){
      sel.k<-info.deviations$Deviations$Cluster==k &
        TypeCase=="Noise"
      sel.core.k<-info.deviations$Deviations$Cluster==k &
        TypeCase=="Core"
      check.cores<-sum(sel.core.k)
      if(check.cores>0){
        dist.k<-Distance[sel.k,sel.core.k,drop=F]<=sel.noise
        forborder_k<-apply(dist.k,1,sum,na.rm=T)
        TypeCase[sel.k][forborder_k>0]<-"Border"
      }}
    cases.values<-TypeCase
  } # closes ref.noise
  
  return(cases.values)
  
}
