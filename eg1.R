rm(list=ls())
library(magrittr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(reshape2)
library(doParallel)
library(parallel)
library(foreach)
source('/Users/zhengxiangyu/Documents/PM2.5/AirQualityAlerts/2_Modeling/Functions/Dynamic_Synthetic_Control_new.R')
setwd('/Users/zhengxiangyu/Documents/PM2.5/AirQualityAlerts/2_Modeling/eg_1')
## 0. read data
{
  data = read.csv('data_eps1.csv')
  data$time = as.POSIXct(data$time)
  # variables settin
  {
    var_x_exo=c('WSPM','humi','dewp','pres')
    var_x_en=c('pm25_lag1')
    var_x=c(var_x_en,var_x_exo)
    var_y=c('pm25')
    var_all = c('id_eps','hour_eps',var_x,var_y)
  }
  # the basic settings 
  {
    N_treated = length(unique(data[data$alert_if==1,'id_eps'])) 
    N_control = length(unique(data$id_eps))-N_treated
    T0=48
    T=T0+24
    data_control=data%>%filter(id_eps>N_treated)%>%dplyr::select(var_all)
    data_treated=data%>%filter(id_eps<=N_treated)%>%dplyr::select(var_all)
    data_treated_mean = aggregate(data_treated[,c(var_x, var_y)],by=list(data_treated$hour_eps),mean)
    colnames(data_treated_mean)[1] = 'hour_eps'
  }
}
## 1. apply the dynamic synthetxic control methods
{
  # 1.1 apply the dynamic synthetic control
  {
    # the modeling
    res_syn=dynamic_syn_con_new(data_control = data_control, data_treated = data_treated_mean, var_y=var_y,
                                var_x_exo = var_x_exo, var_x_en = var_x_en, 
                                T0=T0, T=T, N=N_control)
    ## save the results
    Y0_pre=res_syn$Y0_pre
    Y_observe_mean=data_treated_mean$pm25
    weight_mat = res_syn$weighting_matrix
    ## write results
    {
      res_effect=data.frame(effect = Y_observe_mean-Y0_pre, observe=Y_observe_mean, potential=Y0_pre)%>%mutate(t=1:T)
      mean(res_effect$effect[res_effect$t>48])
      mu1=mean(res_effect$observe[res_effect$t>48])
      mu0=mean(res_effect$potential[res_effect$t>48])
      mu1
      mu0
      mu1-mu0
      1-mu1/mu0
      sqrt(var(res_effect$observe[res_effect$t>48]-res_effect$potential[res_effect$t>48])/24)
      res_effect_summ=data.frame(mu1=mu1,mu0=mu0,diff=mu1-mu0,diff_prop=1-mu1/mu0,T_matched=res_syn$T_matched)
      write.csv(res_effect_summ,'results/effect.csv',row.names = FALSE)
    }
  }
  # 1.3 plot the results for synthetic control estimates
  {
    pre_len=T0
    tr_len=T-T0
    # 2.1 plot the observed average v.s. the synthetic control
    {
      datap=data.frame(hour_eps=1:T, estimated_control=res_syn$Y0_pre,observed_treated=data_treated_mean$pm25)
      datap_long=reshape2::melt(datap, id='hour_eps')
      g1=ggplot(datap_long)+geom_line(aes(x=hour_eps, y=value, color=variable,linetype=variable))+
        ylab(TeX('PM_{2.5} ($\\mu g/m^{3}$)'))+#xlab(TeX('Hours'))+
        theme_bw()+
        scale_linetype_manual(name="Guide1",values= c('longdash', 'solid'),labels = c('Synthetic Control          ','Observed Values'))+
        scale_color_manual(name='Guide1', values=c('brown2','blue'),labels = c('Synthetic Control          ','Observed Values'))+
        annotate("segment", x = 46.2, xend = 48, y = 30, yend = 30, arrow = arrow( length = unit(0.2,"cm")))+
        annotate("text", x = 38, y = 30,label='Air pollution alert',size=6)+
        theme(legend.position = c(0.12,0.92),legend.title = element_blank(),legend.text = element_text(size=10))+
        scale_x_continuous(expand = c(0.01,0.001),breaks = c(1,pre_len/2,pre_len,pre_len+tr_len),labels = c('00:00\n 2016-11-15','00:00\n 2016-11-16','00:00\n 2016-11-17','23:00\n 2016-11-17'))+
        scale_y_continuous(breaks = c(0,100,200,300),limits = c(min(50,min(datap_long$value)),220))+
        geom_vline(xintercept = pre_len, color='dimgrey', linetype="dashed")+
        annotate("text", x= 24, y=180 ,label= 'Pre-intervention',size=6)+
        annotate("text", x= 61, y=180 ,label= 'Post-intervetion',size=6)+
        theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),axis.title.x = element_blank(),axis.title.y = element_text(size=15, hjust=0.5, vjust=1))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        theme(plot.margin = unit(c(0.5,1.5,1,0.5), "cm"))
      ggsave('results/Fig_2_plot_potential_control.pdf',g1,width = 10, height=6)
    }
    # 2.2 plot the gaps between the synthetic control and observed average
    {
      datap=data.frame(hour_eps=1:T, estimated_control=res_syn$Y0_pre,observed_treated=data_treated_mean$pm25)
      datap$value=datap$observed_treated-datap$estimated_control
      g2=ggplot(datap)+geom_line(aes(x=hour_eps, y=value))+
        annotate("segment", x = 46.2, xend = 48, y = -30, yend = -30, arrow = arrow( length = unit(0.3,"cm")))+
        annotate("text", x = 38, y = -30,label='Air pollution alert',size=6)+
        ylab(TeX('PM_{2.5} ($\\mu g/m^{3}$)'))+#xlab(TeX('Hours'))+
        xlab(TeX('Hours'))+
        theme_bw()+
        theme(legend.position = 'top',legend.title = element_blank())+
        scale_x_continuous(expand = c(0.01,0.001),breaks = c(1,pre_len/2,pre_len,pre_len+tr_len),labels = c('00:00\n 2016-11-15','00:00\n 2016-11-16','00:00\n 2016-11-17','23:00\n 2016-11-17'))+
        scale_y_continuous(breaks=c(-60,-40,-20,0,20), limits = c(min(datap$value),40))+
        theme(legend.text = element_text(size=15))+
        geom_vline(xintercept = pre_len, color='dimgrey', linetype="dashed")+
        geom_hline(yintercept = 0, color='dimgrey', linetype="dashed")+
        annotate("text", x= T0/2, y=33 ,label= 'Pre-intervention',size=6)+
        annotate("text", x= T0+12, y=33 ,label= 'Post-intervention',size=6)+
        theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),axis.title.x = element_blank(),axis.title.y = element_text(size=15, hjust=0.65, vjust=1))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        theme(plot.margin = unit(c(0.5,1.5,0.5,0.5), "cm"))
      ggsave('results/Fig_3_plot_estimated_effect.pdf',g2,width = 10, height=6)
    }
  }
}
## 2. assess the identification assumption
{
  # sd of eta_t is proportional to eta_t
  w1=apply(weight_mat^2,2,sum,na.rm=TRUE)[-1]
  v_eta_t = sqrt(1/N_treated + w1)
  # get paramt for rho_t and calculate p-values
  {
    # calculate Bete_{t} and Alpha_{t} and sigma [Y_{it}(0)=beta_{t}X_{it}+alpha_{t}Y_{it}(0)]
    p=length(var_x)
    paramt = data.frame(matrix(0, nrow = T, ncol=p+1))
    colnames(paramt) = c('delt_t',var_x)
    SST=0
    for(t in 1:T0){
      data_temp = data%>%filter(hour_eps==t)%>%dplyr::select(c(var_x, var_y))
      res_ols=lm(pm25~.,data = data_temp)
      paramt[t,]=res_ols$coefficients
      SST=SST+sum(res_ols$residuals^2)
    }
    sigma = sqrt(SST/(nrow(data_temp)*T0))
  }
  # calculate the T0 p-values
  {
    z_diff=((Y_observe_mean-Y0_pre)/(sigma*v_eta_t))[1:T0]
    p_values=2*(1-pnorm(abs(z_diff)))
    plot(sort(p_values))
    # plot p-values
    {
      datap=data.frame(hour_eps=1:T0, pvalues=sort(p_values))
      g=ggplot(datap)+geom_line(aes(x=hour_eps, y=pvalues))+
        geom_point(aes(x=hour_eps, y=pvalues),color='gray')+
        ylab(TeX('p-value'))+#xlab(TeX('Hours'))+
        geom_hline(yintercept = min(p_values), color='dimgrey', linetype="dashed")+
        theme_bw()+
        scale_x_continuous(breaks = c(1,seq(12,48, by=12)))+
        scale_y_continuous(limits=c(0,1),breaks = c(round(min(p_values),3),0.5, 1))+
        theme(axis.title.x = element_blank(),panel.grid.minor = element_blank())
      ggsave('results/Fig_4_pvalues.pdf',g,width = 5.5, height=4)
    }
  }
}
## 3. normalized placebo test
{
  # 3.1 cut data and reserve the 12-hour pretreatment data
  {
    data_plb=data%>%filter(hour_eps>36)
    data_plb$hour_eps=data_plb$hour_eps-36
    data=data_plb
    T0=12
    T=T0+24
    data_control=data.frame(cbind(id_eps=data$id_eps,hour_eps=data$hour_eps,data[,var_x],data[var_y]))%>%filter(id_eps>N_treated)
    colnames(data_control)=c('id_eps','hour_eps',var_x,var_y)
    data_treated=data.frame(cbind(id_eps=data$id_eps,hour_eps=data$hour_eps,data[,var_x],data[var_y]))%>%filter(id_eps<=N_treated)
    colnames(data_treated)=c('id_eps','hour_eps',var_x,var_y)
    data_treated_mean = aggregate(data_treated[,c(var_x, var_y)],by=list(data_treated$hour_eps),mean)
    colnames(data_treated_mean)[1] = 'hour_eps'
  }
  # 3.2 the synthetic control estimates on Beijing
  {
    # the modeling
    res_syn=dynamic_syn_con_new(data_control = data_control, data_treated = data_treated_mean, var_y=var_y,
                                var_x_exo = var_x_exo, var_x_en = var_x_en, 
                                T0=T0, T=T, N=N_control)
    ## save the results
    Y0_pre=res_syn$Y0_pre
    Y_observe_mean=data_treated_mean$pm25
  }
  # 3.3 the placebo trials: run B times
  {
    set.seed(1234)
    parallel_num=detectCores()
    cl=makeCluster(parallel_num)
    registerDoParallel(cl)
    data_plb = data_control
    N_plb = N_control-N_treated
    B=500
    time1=Sys.time()
    bootresult<-foreach(b=1:B,.packages=c('dplyr','LowRankQP','NlcOptim'),.combine=rbind,.inorder=TRUE)%dopar%{
      id_tr_plb = sample((N_treated+1):(N_treated+N_control), N_treated)
      id_co_plb = setdiff((N_treated+1):(N_treated+N_control), id_tr_plb)
      data_treated_plb = data_plb%>%filter(id_eps%in%id_tr_plb)
      data_control_plb = data_plb%>%filter(id_eps%in%id_co_plb)
      N_plb = N_control-N_treated
      
      # summarize data_treated
      {
        mean.na<-function(x){
          return(mean(x, na.rm = TRUE))
        }
        ## variables: 2(id_eps, hour_eps)+p(covariates)+1(y)
        data_treated_mean_plb = aggregate(data_treated_plb[,c(var_x, var_y)],by=list(data_treated_plb$hour_eps),mean.na)
        colnames(data_treated_mean_plb)[1] = 'hour_eps'
      }
      # the modeling
      res_syn_plb=dynamic_syn_con_new(data_control = data_control_plb, data_treated = data_treated_mean_plb, var_y=var_y,
                                      var_x_exo = var_x_exo, var_x_en = var_x_en,
                                      t_begin = 1, T0=T0, T=T,  N=N_plb)
      ## save the results
      bootres_temp=cbind(b,data_treated_mean_plb$pm25, res_syn_plb$Y0_pre, t(data.matrix(res_syn_plb$weighting_matrix[,-1])))
      colnames(bootres_temp)=c('boot_run','observe_mean','synthetic',paste0('unit=',res_syn_plb$weighting_matrix[,1]))
      bootres_temp
    }
    time2=Sys.time()
    time2-time1
    # prepare variables to save the results
    {
      bootresult=data.frame(bootresult)
      Y0_pre_plb = data.frame(matrix(bootresult$synthetic, nrow=T, ncol=B))
      Y_observe_plb = data.frame(matrix(bootresult$observe_mean, nrow=T, ncol=B))
      colnames(Y0_pre_plb)=paste0('plb',1:B)
      colnames(Y_observe_plb)=paste0('plb',1:B)
      weight_plb=array(0, dim = c(N_plb, T, B))
      for(b in 1:B){
        weight_plb[,,b]=t(bootresult[bootresult$boot_run==b,-c(1:3)])
      }
    }
  }
  # 3.4 the normalized placebo test 
  {
    ## obtain the normalized differences
    {
      Y_diff_plb = data.frame(Y0_pre_plb-Y_observe_plb)%>%mutate(hour_eps = 1:T)
      rho = paramt$pm25_lag1
      ## calculate the normalizing coefficients
      {
        p_rescale<-function(weights_mat, t, rho, T0, N_treated){
          pt=0
          if(t<=T0+1){
            pt=(1/N_treated+sum(weights_mat[,t]^2, na.rm = TRUE))
          }else if(t>=T0+2){
            pt=(1/N_treated+sum(weights_mat[,t]^2, na.rm = TRUE))
            for(s in (T0+1):(t-1)){
              pt=pt+prod(rho[(s+1):t])^2*(1/N_treated+sum(weights_mat[,s]^2))
            }
          }
          return(sqrt(pt))
        }
        p_matrix = array(0, dim = c(T, B))
        for(b in 1:B){
          weights_mat=weight_plb[,,b]
          for(t in 1:T){
            if(t<=T0){
              #p_matrix[b,t]=1
              p_matrix[t, b]=p_rescale(weights_mat, t, rho, T0, N_treated)
            }else{
              p_matrix[t, b]=p_rescale(weights_mat, t, rho, T0, N_treated)
            }
          }
        }
        p_treated = rep(1,T)
        weights_mat=res_syn$weighting_matrix[,-1]
        for(t in 1:T){
          if(t<=T0){
            p_treated[t]=p_rescale(weights_mat, t, rho, T0, N_treated)
          }else{
            p_treated[t]=p_rescale(weights_mat, t, rho, T0, N_treated)
          }
        }
        p_coef=(p_treated%*%t(rep(1, B))/p_matrix)
      }
      ## normalize the differences
      {
        Y_diff_plb_norm = Y_diff_plb
        Y_diff_plb_norm[,1:B] = Y_diff_plb_norm[,1:B]*p_coef
      }
    }
    ## plot the normalized results 
    {
      pre_len=T0
      tr_len=T-T0
      datap=Y_diff_plb_norm
      datap=datap[,intersect(names(which(apply(datap,2,min)>(-100))),names(which(apply(datap,2,max)<(100))))]
      datap_long
      datap_long=reshape2::melt(datap, id='hour_eps')
      datap_long$type='placebo'
      datap_add = data.frame(hour_eps=1:T, value=Y0_pre-Y_observe_mean, 
                             variable = 'treated', type='treated')
      datap_long=rbind(datap_long, datap_add)
      datap_long$type=factor(datap_long$type,levels = c('treated', 'placebo'), labels = c('treated', 'placebo'))
      pdf('results/Fig_7_placebo_summary_norm.pdf', width=9,heigh=6.5)
      g_summ=ggplot(data=datap_long)+
        theme_bw()+
        geom_line(data=datap_long[datap_long$type=='placebo',],aes(x=hour_eps, y=-value,group=variable),color='gray85',alpha=0.5)+
        geom_line(data=datap_long[datap_long$type=='treated',],aes(x=hour_eps, y=-value,group=variable),color='black')+
        geom_vline(xintercept = pre_len, color='dimgrey', linetype="dashed")+
        ylab(TeX('PM_{2.5} ($\\mu g/m^{3}$)'))+#xlab(TeX('Hours'))+
        annotate("segment", x = 11, xend = 12, y = -60, yend = -60, arrow = arrow( length = unit(0.2,"cm")))+
        annotate("text", x = 7, y = -60,label='Air pollution alert',size=5)+
        scale_x_continuous(expand = c(0,0),breaks = c(1,pre_len,pre_len+tr_len),labels = c('12:00\n 2016-11-16','00:00\n 2016-11-17','23:00\n 2016-11-17   '))+
        scale_y_continuous(limits=c(-100,100),breaks = c(-80,-40,0,40))+
        geom_line(aes(x=hour_eps, y=-value-1000, color=type, group=variable))+
        scale_color_manual(values=c('black','gray85'),labels=c('Beijing','Placebo'))+
        theme(legend.position = c(0.07,0.93),legend.title = element_blank(),legend.text = element_text(size=10))+
        theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=18),axis.title.x = element_blank(),axis.title.y = element_text(size=15, hjust=0.5, vjust=1))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        theme(plot.margin = unit(c(0.5,1.5,1.2,0.5), "cm"))
      print(g_summ)
      dev.off()
    }
    ## calculate the pvalues
    {
      Y_diff_treated=Y0_pre-Y_observe_mean
      pvalues=apply(Y_diff_plb_norm[,1:B]>(Y_diff_treated)%*%rep(1,B),1,mean)[(T0+1):T]
      res_pvalue=data.frame(time=paste0('$T_{0}+',1:(T-T0),'$'), pvalues=pvalues)
      write.csv(t(res_pvalue), 'results/pvalues_placebo_random.csv')
      res=t((read.csv('/Users/zhengxiangyu/Documents/PM2.5/AirQualityAlerts/2_Modeling/eg_1/results/pvalues_placebo_random.csv'))[,-1])
      pvalues=data.frame(hour_eps=1:24, pvalue=as.numeric(res[,2]))
      pvalues$rank=rank(pvalues$pvalue)
      pvalues$level= round(pvalues$rank/24*0.05,2)
      pvalues$significant=pvalues$pvalue<= pvalues$level
    }
  }
}
