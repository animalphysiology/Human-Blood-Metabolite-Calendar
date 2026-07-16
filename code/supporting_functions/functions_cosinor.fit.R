
cosinor.model.test<-function(m,features,metadata) {
  y=features[m,metadata$Sample]|>log2()|>unlist()
  suppressMessages(
    { fit.cosinor<-lmer(y~sin.a+cos.a+Sex+(1|Subject),data=metadata,REML=F)}
  )
  para<-fixef(fit.cosinor)
  return(para)
}

cosinor.sd.model.test<-function(m,features,metadata) {
  y=features[m,metadata$Sample]|>log2()|>unlist()
  fit.cosinor.sd<-lmer(y~(sin.a+cos.a)*Sex+(1|Subject),
                       data=metadata,REML=F)
  para<-fixef(fit.cosinor.sd)
  return(para)
}




# bootstrap parameter calculation by shuffling subjects
# from lmer object
subjects.boot.cosinor<-function(fit.lmm,metadata,n=100,seed=777){
  subjects = unique(fit.lmm@frame$Subject)
  m=nrow(metadata)
  set.seed(seed)
  lapply(1:n,function(x){
    df.boot=left_join(data.frame("Subject"=sample(subjects,replace = T)),
                      metadata,by=join_by(Subject))
    fit.boot=lmer(formula(fit.lmm),data=df.boot[-sample(1:m,1),])
    return(fixef(fit.boot))
  })
  
}

