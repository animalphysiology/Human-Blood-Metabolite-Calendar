library(data.table)



mz.rt.matching<-function(dt1,dt2,mz.tol=10,rt.tol=NULL,high.mass=0.001){
  dt1 <- as.data.table(human_amide_BGI)
  dt2 <- as.data.table(human_amide_BGI_2nd)

  dt1[, `:=`(
    l = CalcMz - mz.tol * CalcMz / 1e6,
    u = CalcMz + mz.tol * CalcMz / 1e6
  )]
  
  dt2[, `:=`(
    l = CalcMz ,
    u = CalcMz 
  )]
  
  setkey(dt1, Ion.Mode, l, u)
  setkey(dt2, Ion.Mode, l, u)
  
  match.df <- foverlaps(
    dt2, dt1,
    by.x = c("Ion.Mode", "l", "u"),
    by.y = c("Ion.Mode", "l", "u"),
    type = "within",
    nomatch = 0L
  )
  
  match.df[, `:=`(
    d.ppm = (CalcMz - i.CalcMz) / i.CalcMz * 1e6,
    d.rt  = i.RT - RT
  )]
  
  match.df <- match.df[, .(
    ref.id   = unique_id,
    id         = i.unique_id,
    mode       = Ion.Mode,
    name       = i.Name,
    ref.name = Name,
    mz         = i.CalcMz,
    ref.mz   = CalcMz,
    d.ppm,
    ref.rt   = RT,
    rt         = i.RT,
    d.rt
  )]
  
  match.df.all<-match.df[(abs(match.df$d.ppm)<mz.tol|abs(match.df$d.ppm)<=high.mass),]
  
  if(!is.null(rt.tol)){
    match.df.all<-match.df.all[(abs(match.df.all$d.rt)<rt.tol),]
  }
  
  return(match.df.all)
}



timetable<- function(SOMs,metadata,features){
  timetable.pred<-metadata
  rn=intersect(SOMs$unique_id,rownames(r))
  SOMs<-SOMs[rn,]
  for (month in metadata$Sample) {
    y=features[SOMs$unique_id,month]|>unlist()
    phase=SOMs[,paste0("phase")]
    fit.lm<-glm(y~sin(phase*w)+cos(phase*w))
    para<-coef(fit.lm)
    para["RS"]=1-deviance(fit.lm)/deviance(glm(y~1))
    para["amplitude"]<-2*sqrt(para["sin(phase * w)"]^2+para["cos(phase * w)"]^2)
    para["predicted"]=
      (atan2(para["sin(phase * w)"],para["cos(phase * w)"])/w)%%period
    para["actual"]=metadata[metadata$Sample==month,]$day_of_year%%period
    para["error"]=circular.diff(para["predicted"],para["actual"])
    
    timetable.pred[metadata$Sample==month,names(para)] <- para
  }
  return(timetable.pred)
}