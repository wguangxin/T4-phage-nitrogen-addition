library(readr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(patchwork)


stat.df <- read_csv("data/stat.csv", show_col_types=FALSE) %>%
  mutate(
    N_f=factor(N_add, levels=c(0,2,5,10,20,50),
               labels=c("N0","N2","N5","N10","N20","N50")),
    Year=factor(Year), Block=factor(Block), Plot_ID=factor(Plot)
  )



p_format <- function(p) if (p < 0.001) "< 0.001" else paste0("= ",sprintf("%.3f",p))

make_panel <- function(data, variable, title, ylab, log10_transform=FALSE) {
  model <- lmer(.response ~ N_f*Year + (1|Block), data=data, REML=TRUE)
  test <- anova(model, type=3, ddf="Kenward-Roger")
  emm <- emmeans(model, ~ N_f | Year)
  letters.df <- multcomp::cld(emm, adjust="tukey", Letters=letters,
                              sort=FALSE, reversed=TRUE) %>%
    as.data.frame() %>%
    mutate(.group=gsub(" ","",.group))
  offset <- 0.055*diff(range(data$.response))
  letter.pos <- data %>% group_by(Year,N_f) %>%
    summarise(label_y=max(.response)+offset, .groups="drop")
  letters.df <- left_join(letters.df, letter.pos, by=c("Year","N_f"))
  label <- paste0("N: P ",p_format(test["N_f","Pr(>F)"]),
                  "\nYear: P ",p_format(test["Year","Pr(>F)"]),
                  "\nN × Year: P ",p_format(test["N_f:Year","Pr(>F)"]))
  p.df <- data.frame(Year=factor("2016",levels=levels(data$Year)),
                     x=1, y=Inf, label=label)
  
  p <- ggplot(data, aes(N_f,.response)) +
    annotate("rect",xmin=5.5,xmax=6.5,ymin=-Inf,ymax=Inf,
             fill="#F4DFA3",alpha=0.35) +
    geom_point(aes(color=Year),position=position_jitter(width=0.07,height=0),
               size=2,alpha=0.65,show.legend=FALSE) +
    stat_summary(fun.data=mean_se,geom="errorbar",width=0.15,linewidth=0.55) +
    stat_summary(fun=mean,geom="point",shape=23,size=2,
                 fill="white",color="black",stroke=0.85) +
    geom_text(data=letters.df,aes(x=N_f,y=label_y,label=.group,color=Year),
              inherit.aes=FALSE,size=3.5,fontface="bold",show.legend=FALSE) +
    geom_text(data=p.df,aes(x=x,y=y,label=label),inherit.aes=FALSE,
              hjust=0,vjust=1.15,size=2.8,lineheight=1.05) +
    facet_wrap(~Year,nrow=1) +
    scale_color_manual(values=c("2016"="#3C78A8","2018"="#D95F59")) +
    scale_y_continuous(expand=expansion(mult=c(0.04,0.18))) +
    labs(title=title,
         x=expression("Nitrogen addition (g N m"^{-2}*" yr"^{-1}*")"),
         y=ylab) +
    theme_classic(base_size=9) +
    theme(
      plot.title=element_text(size=10,face="bold"),
      strip.background=element_rect(fill="grey93",color="grey35",linewidth=0.4),
      strip.text=element_text(face="bold"),
      axis.text=element_text(color="black"),
      axis.title=element_text(color="black"),
      panel.spacing=grid::unit(4,"mm")
    )
  
  list(plot=p,model=model,test=test,emmeans=emm,letters=letters.df)
}

g20 <- make_panel(stat.df,"g20.rab","A  g20",
                  expression(log[10]*" g20 relative abundance"),TRUE)
g23 <- make_panel(stat.df,"g23.rab","B  g23",
                  expression(log[10]*" g23 relative abundance"),TRUE)
host <- make_panel(stat.df,"host.rab","C  Potential host",
                   expression(log[10]*" potential host relative abundance"),TRUE)
richness <- make_panel(stat.df,"host_Shannon","D  Potential host",
                       "Potential host Shannon index",FALSE)




library(readr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)

stat.df <- read_csv("data/stat.csv",show_col_types=FALSE) %>%
  mutate(N_f=factor(N_add,levels=c(0,2,5,10,20,50),
                    labels=c("N0","N2","N5","N10","N20","N50")),
         Year=factor(Year),Block=factor(Block))

fig.dir <- "figure/"

p_format <- function(p) if(p<0.001) "< 0.001" else paste0("= ",sprintf("%.3f",p))

make_bar_plot <- function(data,variable,y_label,log_value=TRUE){
  
  model <- lmer(.response ~ N_f*Year + (1|Block),data=data,REML=TRUE)
  test <- anova(model,type=3,ddf="Kenward-Roger")
  emm <- emmeans(model,~N_f|Year)

  bar.df <- data %>% group_by(Year,N_f) %>%
    summarise(mean=mean(.response),se=sd(.response)/sqrt(n()),.groups="drop") %>%
    mutate(x=as.numeric(N_f))

  letter.df <- multcomp::cld(emm,adjust="tukey",Letters=letters,
                             sort=FALSE,reversed=TRUE) %>%
    as.data.frame() %>% mutate(.group=gsub(" ","",.group))
  offset <- 0.055*diff(range(data$.response))
  letter.pos <- data %>% group_by(Year,N_f) %>%
    summarise(label_y=max(.response)+offset,.groups="drop")
  letter.df <- left_join(letter.df,letter.pos,by=c("Year","N_f")) %>%
    mutate(x=as.numeric(N_f))
  
  label <- paste0("N: P ",p_format(test["N_f","Pr(>F)"]),
                  "\nYear: P ",p_format(test["Year","Pr(>F)"]),
                  "\nN × Year: P ",p_format(test["N_f:Year","Pr(>F)"]))
  label.df <- data.frame(Year=factor("2016",levels=levels(data$Year)),
                         x=0.68,y=Inf,label=label)
  
  baseline <- if(log_value) floor(min(data$.response)) else 0
  
  p <- ggplot() +
    annotate("rect",xmin=5.5,xmax=6.5,ymin=-Inf,ymax=Inf,
             fill="#F4DFA3",alpha=0.35) +
    geom_rect(data=bar.df,aes(xmin=x-0.33,xmax=x+0.33,
                              ymin=baseline,ymax=mean,fill=Year),
              color="black",linewidth=0.3,alpha=0.8,show.legend=FALSE) +
    geom_errorbar(data=bar.df,aes(x=x,ymin=mean-se,ymax=mean+se),
                  width=0.16,linewidth=0.55) +
    #geom_point(data=data,aes(x=as.numeric(N_f),y=.response),
    #           position=position_jitter(width=0.08,height=0),
    #           shape=21,fill="white",color="black",size=1.7,stroke=0.4) +
    geom_text(data=letter.df,aes(x=x,y=label_y,label=.group,color=Year),
              size=3.5,fontface="bold",show.legend=FALSE) +
    geom_text(data=label.df,aes(x=x,y=y,label=label),
              hjust=0,vjust=1.15,size=2.8,lineheight=1.05) +
    facet_wrap(~Year,nrow=1) +
    scale_x_continuous(breaks=1:6,labels=levels(data$N_f),limits=c(0.5,6.5)) +
    scale_y_continuous(limits=c(baseline,NA),expand=expansion(mult=c(0,0.22))) +
    scale_fill_manual(values=c("2016"="#3C78A8","2018"="#D95F59")) +
    scale_color_manual(values=c("2016"="#3C78A8","2018"="#D95F59")) +
    labs(x=expression("Nitrogen addition (g N m"^{-2}*" yr"^{-1}*")"),
         y=y_label) +
    theme_classic(base_size=9) +
    theme(strip.background=element_rect(fill="grey93",color="grey35",linewidth=0.4),
          strip.text=element_text(face="bold"),
          axis.text=element_text(color="black"),
          axis.title=element_text(color="black"),
          panel.spacing=grid::unit(5,"mm"))
  
  list(plot=p,model=model,test=test,letters=letter.df,bar_data=bar.df)
}

# 分别生成四张图
g20 <- make_bar_plot(stat.df,"g20.rab",
                     expression(log[10]*" g20 relative abundance"))
g23 <- make_bar_plot(stat.df,"g23.rab",
                     expression(log[10]*" g23 relative abundance"))
host <- make_bar_plot(stat.df,"host.rab",
                      expression(log[10]*" host relative abundance"))
richness <- make_bar_plot(stat.df,"host_Shannon",
                          "Host Shannon index",log_value=FALSE)





library(readr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(patchwork)

stat.df <- read_csv("data/stat.csv",show_col_types=FALSE) %>%
  mutate(
    N_f=factor(N_add,levels=c(0,2,5,10,20,50),
               labels=c("N0","N2","N5","N10","N20","N50")),
    Year=factor(Year),Block=factor(Block),
    VHR_g20=g20.rab/host.rab,
    VHR_g23=g23.rab/host.rab,
    VHR_g20_log=log10(VHR_g20),
    VHR_g23_log=log10(VHR_g23)
  )


p_format <- function(p) if(p<0.001) "< 0.001" else paste0("= ",sprintf("%.3f",p))

plot_vhr <- function(data,response,gene){
  data$.response <- data[[response]]
  model <- lmer(.response ~ N_f*Year + (1|Block),data=data,REML=TRUE)
  test <- anova(model,type=3,ddf="Kenward-Roger")
  emm <- emmeans(model,~N_f|Year)
  letter.df <- multcomp::cld(emm,adjust="tukey",Letters=letters,
                             sort=FALSE,reversed=TRUE) %>%
    as.data.frame() %>% mutate(.group=gsub(" ","",.group))
  
  offset <- 0.055*diff(range(data$.response))
  letter.pos <- data %>% group_by(Year,N_f) %>%
    summarise(label_y=max(.response)+offset,.groups="drop")
  letter.df <- left_join(letter.df,letter.pos,by=c("Year","N_f"))
  mean.df <- data %>% group_by(Year,N_f) %>%
    summarise(mean=mean(.response),se=sd(.response)/sqrt(n()),.groups="drop")
  
  label <- paste0("N: P ",p_format(test["N_f","Pr(>F)"]),
                  "\nYear: P ",p_format(test["Year","Pr(>F)"]),
                  "\nN × Year: P ",p_format(test["N_f:Year","Pr(>F)"]))
  label.df <- data.frame(Year=factor("2016",levels=levels(data$Year)),
                         x=1,y=Inf,label=label)
  
  p <- ggplot(data,aes(N_f,.response,color=Year)) +
    geom_point(position=position_jitter(width=0.055,height=0),
               size=2,alpha=0.45,show.legend=FALSE) +
    geom_line(data=mean.df,aes(x=N_f,y=mean,color=Year,group=Year),
              linewidth=0.85,inherit.aes=FALSE) +
    geom_errorbar(data=mean.df,aes(x=N_f,ymin=mean-se,ymax=mean+se,color=Year),
                  width=0.13,linewidth=0.6,inherit.aes=FALSE) +
    geom_point(data=mean.df,aes(x=N_f,y=mean,color=Year),
               shape=23,fill="white",size=3.3,stroke=0.9,inherit.aes=FALSE) +
    geom_text(data=letter.df,aes(x=N_f,y=label_y,label=.group,color=Year),
              inherit.aes=FALSE,size=3.4,fontface="bold",show.legend=FALSE) +
    geom_text(data=label.df,aes(x=x,y=y,label=label),
              inherit.aes=FALSE,hjust=0,vjust=1.15,size=2.8,lineheight=1.05) +
    facet_wrap(~Year,nrow=1) +
    scale_color_manual(values=c("2016"="#3C78A8","2018"="#D95F59")) +
    scale_y_continuous(expand=expansion(mult=c(0.06,0.20))) +
    labs(x=expression("Nitrogen addition (g N m"^{-2}*" yr"^{-1}*")"),
         y=bquote(log[10]~.(gene)~"/potential-host ratio")) +
    theme_classic(base_size=9) +
    theme(legend.position="none",
          strip.background=element_rect(fill="grey93",color="grey35",linewidth=0.4),
          strip.text=element_text(face="bold"),
          axis.text=element_text(color="black"),
          axis.title=element_text(color="black"),
          panel.spacing=grid::unit(5,"mm"))
  
  list(plot=p,model=model,test=test,letters=letter.df,means=mean.df)
}

g20.vhr <- plot_vhr(stat.df,"VHR_g20_log","g20")
g23.vhr <- plot_vhr(stat.df,"VHR_g23_log","g23")




library(readr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)

stat.df <- read_csv("data/stat.csv",show_col_types=FALSE) %>%
  mutate(N_f=factor(N_add,levels=c(0,2,5,10,20,50),
                    labels=c("N0","N2","N5","N10","N20","N50")),
         Year=factor(Year),Block=factor(Block))

fig.dir <- "figure"

variables <- c(g20="g20.rab",g23="g23.rab",
               Proteobacteria="Proteobacteria.rab",
               Cyanobacteria="Cyanobacteria.rab",
               `Overall potential host`="host.rab")

get_effect <- function(data,variable,label){
  data$.response <- log10(data[[variable]])
  
  # 与前面分析一致：Block为随机效应
  model <- lmer(.response ~ N_f*Year + (1|Block),data=data,REML=TRUE)
  emm <- emmeans(model,~N_f|Year)
  con <- contrast(emm,list("N50 vs N0–N20"=c(-0.2,-0.2,-0.2,-0.2,-0.2,1)))
  out <- as.data.frame(summary(con,infer=TRUE,adjust="none"))
  out %>% transmute(Indicator=label,Year,
                    effect=estimate/log10(2),
                    lower=lower.CL/log10(2),
                    upper=upper.CL/log10(2),
                    p.value=p.value)
}

effect.df <- bind_rows(lapply(names(variables),function(nm)
  get_effect(stat.df,variables[[nm]],nm)))
indicator.order <- c("Overall potential host","Cyanobacteria",
                     "Proteobacteria","g23","g20")
effect.df <- effect.df %>%
  mutate(Indicator=factor(Indicator,levels=indicator.order),
         Year=factor(Year,levels=c("2016","2018")),
         ypos=as.numeric(Indicator)+ifelse(Year=="2016",-0.13,0.13),
         cap_low=ypos-0.045,cap_high=ypos+0.045)

print(effect.df %>% dplyr::select(Indicator,Year,effect,lower,upper,p.value))

p.effect <- ggplot(effect.df,aes(color=Year)) +
  geom_vline(xintercept=0,linetype=2,color="grey50",linewidth=0.5) +
  geom_segment(aes(x=lower,xend=upper,y=ypos,yend=ypos),linewidth=0.9) +
  geom_segment(aes(x=lower,xend=lower,y=cap_low,yend=cap_high),linewidth=0.9) +
  geom_segment(aes(x=upper,xend=upper,y=cap_low,yend=cap_high),linewidth=0.9) +
  geom_point(aes(x=effect,y=ypos,shape=Year),size=3.4,stroke=0.8) +
  scale_shape_manual(values=c("2016"=16,"2018"=17)) +
  scale_color_manual(values=c("2016"="#3C78A8","2018"="#D95F59")) +
  scale_y_continuous(breaks=1:5,labels=indicator.order,
                     limits=c(0.5,5.5),expand=c(0,0)) +
  scale_x_continuous(expand=expansion(mult=c(0.08,0.08))) +
  labs(x=expression(log[2]*" fold change: N50 vs N0–N20"),y=NULL,
       color=NULL,shape=NULL) +
  theme_bw(base_size=10) +
  theme(
    panel.border=element_rect(color="black",fill=NA,linewidth=0.6),
    panel.grid=element_blank(),
    axis.text=element_text(color="black"),
    axis.title=element_text(color="black"),
    legend.position="inside",
    legend.position.inside=c(0.88,0.20),
    legend.background=element_blank(),
    plot.margin=margin(6,10,6,6)
  )

print(p.effect)




library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)

#-------------------- 1. 新数据与变量 --------------------
stat.df <- read_csv("data/stat.csv",show_col_types=FALSE)

required <- c("Block","Year","g20.rab","g23.rab","host.rab","Ph","AN",
              "NO3_N","NH4_N","available_P","available_C","total_C",
              "total_N","water_content")

driver.labels <- c(
  Ph="pH",AN_t="AN",NO3_t="NO3-N",NH4_t="NH4-N",
  available_P="AP",available_C="AC",total_C="TC",total_N="TN",
  water_content="Moisture",C_N="C/N"
)
response.labels <- c(VHR_g20_log="g20-VHR",VHR_g23_log="g23-VHR")
pearson.adjust <- "BH"

stat.df <- stat.df %>%
  mutate(
    Block=factor(Block),Year=factor(Year),
    VHR_g20_log=log10(g20.rab/host.rab),
    VHR_g23_log=log10(g23.rab/host.rab),
    C_N=total_C/total_N,
    AN_t=log10(AN+1),NO3_t=log10(NO3_N+1),NH4_t=log10(NH4_N+1)
  ) %>%
  filter(Year %in% c("2016","2018"))

model.control <- lme4::lmerControl(
  check.conv.singular=list(action="ignore",tol=1e-4)
)

analyse_driver <- function(data,period){
  vars <- c(names(response.labels),names(driver.labels))
  d <- data %>%
    filter(if_all(all_of(vars),~is.finite(.x)),!is.na(Block)) %>%
    droplevels() %>% as.data.frame()
  
  bad <- vars[vapply(d[vars],function(z) !is.finite(sd(z)) || sd(z)==0,logical(1))]

  d[vars] <- lapply(d[vars],function(z) as.numeric(scale(z)))
  
  res <- bind_rows(lapply(names(response.labels),function(y){
    bind_rows(lapply(names(driver.labels),function(x){
      model <- lmerTest::lmer(as.formula(paste(y,"~",x,"+ (1|Block)")),
                              data=d,REML=TRUE,control=model.control)
      test <- anova(model,type=3,ddf="Kenward-Roger")
      r2 <- MuMIn::r.squaredGLMM(model)
      data.frame(Period=period,Response=unname(response.labels[y]),
                 Predictor=x,Beta=unname(lme4::fixef(model)[x]),
                 R2m=unname(r2[1,"R2m"]),P=test[x,"Pr(>F)"],N=nrow(d))
    }))
  })) %>%
    group_by(Response) %>%
    mutate(P_FDR=p.adjust(P,method="BH"),
           Strength=sqrt(pmax(R2m,0))) %>%
    ungroup()
  R <- cor(d[names(driver.labels)],method="pearson")
  ij <- which(upper.tri(R),arr.ind=TRUE)
  mat <- data.frame(i=ij[,1],j=ij[,2],r=R[ij])
  mat$P <- vapply(seq_len(nrow(mat)),function(k){
    cor.test(d[[names(driver.labels)[mat$i[k]]]],
             d[[names(driver.labels)[mat$j[k]]]],method="pearson")$p.value
  },numeric(1))
  mat$P_adjusted <- p.adjust(mat$P,method=pearson.adjust)
  mat$Stars <- ifelse(mat$P_adjusted<0.001,"***",
                      ifelse(mat$P_adjusted<0.01,"**",
                             ifelse(mat$P_adjusted<0.05,"*","")))
  
  list(result=res,correlation=R,matrix=mat,period=period)
}

result.all <- analyse_driver(stat.df,"2016 and 2018")
result.2016 <- analyse_driver(filter(stat.df,Year=="2016"),"2016")
result.2018 <- analyse_driver(filter(stat.df,Year=="2018"),"2018")

all.results <- bind_rows(result.all$result,result.2016$result,result.2018$result)

r.max <- max(0.1,ceiling(max(all.results$Strength)*10)/10)

plot_driver <- function(object,r.max){
  curve.colors <- c(Positive="#79AFCB",Negative="#E5AA82")
  thin <- 0.25                        
  significant.width <- c(0.45,1.25)   
  label.size <- 2.8
  number.size <- 2.4
  star.size <- 3.0
  
  R <- object$correlation
  n <- ncol(R)
  vars <- colnames(R)
  labels <- unname(driver.labels[vars])
  
  nodes <- data.frame(Response=c("g20-VHR","g23-VHR"),
                      x=c(-2.0,0.0),y=c(5.0,2.0))
  curve.bend <- 0.10                 
  label.t <- c("g20-VHR"=0.60,"g23-VHR"=0.78)
  number.offset <- 0.28             
  
  mat <- object$matrix
  mat$x <- mat$j
  mat$y <- n-mat$i+1
  mat$half <- 0.43*sqrt(abs(mat$r))
  diagonal <- data.frame(x=seq_len(n)+0.35,y=n-seq_len(n)+1,label=labels)
  top <- data.frame(x=2:n,y=n+0.65,label=labels[-1])
  
  name.gap <- c(Ph=0.65,AN_t=0.65,NO3_t=1.15,NH4_t=1.15,
                available_P=0.65,available_C=0.65,total_C=0.65,
                total_N=0.65,water_content=1.55,C_N=0.85)
  
  links <- object$result %>%
    mutate(i=match(Predictor,vars),
           Direction=ifelse(Beta>=0,"Positive","Negative"),
           Sig=ifelse(P_FDR<0.05,"FDR < 0.05","FDR >= 0.05"),
           Width=ifelse(P_FDR<0.05,
                        significant.width[1]+diff(significant.width)*Strength/r.max,thin)) %>%
    left_join(nodes,by="Response") %>%
    mutate(xend=i+0.35-unname(name.gap[Predictor])-0.15,
           yend=n-i+1,edge=row_number(),
           cx=(x+xend)/2-curve.bend*(yend-y),
           cy=(y+yend)/2+curve.bend*(xend-x))
  paths <- bind_rows(lapply(seq_len(nrow(links)),function(k){
    a <- links[k,]
    t <- seq(0,1,length.out=100)
    data.frame(x=(1-t)^2*a$x+2*(1-t)*t*a$cx+t^2*a$xend,
               y=(1-t)^2*a$y+2*(1-t)*t*a$cy+t^2*a$yend,
               edge=a$edge,Direction=a$Direction,Sig=a$Sig,Width=a$Width)
  }))
  
  line.labels <- links %>%
    filter(P_FDR<0.05) %>%
    mutate(t=unname(label.t[Response]),
           lx=(1-t)^2*x+2*(1-t)*t*cx+t^2*xend,
           ly=(1-t)^2*y+2*(1-t)*t*cy+t^2*yend,
           dx=2*(1-t)*(cx-x)+2*t*(xend-cx),
           dy=2*(1-t)*(cy-y)+2*t*(yend-cy),
           tangent.length=pmax(sqrt(dx^2+dy^2),1e-8),
           lx=lx-number.offset*dy/tangent.length,
           ly=ly+number.offset*dx/tangent.length,
           Label=sprintf("%.2f",Strength))
  
  bottom <- min(mat$y)-0.5
  upper <- max(mat$y)+0.5
  bar.x <- n+0.85               
  edges <- seq(bottom,upper,length.out=301)
  bar <- data.frame(ymin=head(edges,-1),ymax=tail(edges,-1),
                    r=seq(-1,1,length.out=300))
  ticks <- data.frame(r=c(-1,-0.5,0,0.5,1))
  ticks$y <- bottom+(ticks$r+1)/2*(upper-bottom)
  
  r.breaks <- pretty(c(0,r.max),n=3)
  r.breaks <- r.breaks[r.breaks>0 & r.breaks<=r.max]
  width.breaks <- significant.width[1]+diff(significant.width)*r.breaks/r.max
  
  ggplot() +
    geom_tile(data=mat,aes(x,y),width=1,height=1,
              fill="white",color="grey70",linewidth=0.3) +
    geom_rect(data=mat,aes(xmin=x-half,xmax=x+half,ymin=y-half,ymax=y+half,fill=r)) +
    geom_text(data=filter(mat,abs(r)<0.65),aes(x,y,label=Stars),
              size=star.size,color="black",fontface="bold") +
    geom_text(data=filter(mat,abs(r)>=0.65),aes(x,y,label=Stars),
              size=star.size,color="white",fontface="bold") +
    geom_path(data=filter(paths,Sig=="FDR >= 0.05"),
              aes(x,y,group=edge,color=Direction,linewidth=Width,linetype=Sig),
              lineend="round") +
    geom_path(data=filter(paths,Sig=="FDR < 0.05"),
              aes(x,y,group=edge,color=Direction,linewidth=Width,linetype=Sig),
              lineend="round") +
    geom_text(data=line.labels,aes(lx,ly,label=Label,color=Direction),
              size=number.size,show.legend=FALSE) +
    geom_text(data=diagonal,aes(x,y,label=label),hjust=1,size=label.size) +
    geom_text(data=top,aes(x,y,label=label),angle=45,hjust=0,size=label.size) +
    geom_point(data=nodes,aes(x,y),shape=16,size=3,color="grey25") +
    geom_text(data=nodes,aes(x,y,label=Response),nudge_x=-0.25,hjust=1,
              size=3.2,fontface="bold") +
    geom_rect(data=bar,aes(ymin=ymin,ymax=ymax,fill=r),
              xmin=bar.x,xmax=bar.x+0.35) +
    annotate("rect",xmin=bar.x,xmax=bar.x+0.35,ymin=bottom,ymax=upper,
             fill=NA,color="grey30",linewidth=0.3) +
    geom_segment(data=ticks,aes(y=y,yend=y),
                 x=bar.x+0.35,xend=bar.x+0.47,linewidth=0.3) +
    geom_text(data=ticks,aes(y=y,label=r),x=bar.x+0.58,hjust=0,size=2.8) +
    annotate("text",x=bar.x+1.65,y=(bottom+upper)/2,
             label="Pearson's r",angle=90,size=3) +
    
    scale_fill_gradient2(low="#B2182B",mid="white",high="#2166AC",
                         limits=c(-1,1),guide="none") +
    scale_color_manual(values=curve.colors,limits=c("Positive","Negative"),
                       name="LMM association") +
    scale_linetype_manual(
      values=c("FDR < 0.05"="solid","FDR >= 0.05"="dashed"),
      limits=c("FDR < 0.05","FDR >= 0.05"),drop=FALSE,
      name="Significance"
    ) +
    scale_linewidth_identity(
      breaks=width.breaks,labels=signif(r.breaks,2),
      limits=c(thin,significant.width[2]),guide="legend",
      name=expression(sqrt(R[m]^2))
    ) +
    coord_fixed(xlim=c(min(nodes$x)-2.5,n+3.2),
                ylim=c(0,n+4.5),clip="off") +
    labs(title=object$period,
         caption=paste0("Pearson tests (",pearson.adjust,
                        "): * < 0.05; ** < 0.01; *** < 0.001")) +
    theme_void(base_size=9) +
    theme(plot.title=element_text(size=11,face="bold",hjust=0),
          plot.caption=element_text(size=8,hjust=1),
          legend.position="bottom",legend.box="horizontal",
          legend.title=element_text(size=8),legend.text=element_text(size=8),
          legend.key.width=grid::unit(9,"mm"),
          plot.margin=margin(5,10,5,10)) +
    guides(
      color=guide_legend(order=1,override.aes=list(linewidth=0.8,linetype="solid")),
      linewidth=guide_legend(order=2,override.aes=list(color="grey35",linetype="solid")),
      linetype=guide_legend(order=3,override.aes=list(color="grey35",
                                                      linewidth=c(0.8,thin)))
    )
}

