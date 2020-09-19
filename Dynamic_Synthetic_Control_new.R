# 【input】
## data_control (T x N_treated x (1+p)) data_treated (T x N_control x (1+p)) var_x(1 x p) var_y 
# 【output】
## Y_pre (T x 1) weight (T x N_control) approximate (Txp) app.diff (Txp)
require(LowRankQP)
require(NlcOptim)
## t_begin<=T0+1
## the rank of variables in Z1 : var_en var_exo var_loadings
dynamic_syn_con_new<-function(data_control, data_treated, var_y, var_x_exo, var_x_en, T0 , T,  N, N0=1,t_begin=1,K_lag=1, weight_sum=NULL, bound=1){
  var_x_en=var_x_en
  id_eps_ctl=unique(data_control$id_eps)
  id_eps_trt=unique(data_treated$id_eps)
  weight_mat = cbind(data.frame(id_eps_control = unique(data_control$id_eps)), data.frame(matrix(nrow=N,ncol=T)))
  colnames(weight_mat)[2:ncol(weight_mat)] = paste0('t=',1:T)
  var_diff_mat = data.frame()
  var_Z_mat = data.frame()
  var_y_pre = array(0, c(T,1))
  # prepare the V matrix for minimizing the distance
  # which assigns more importance to the variables with higher correlations
  {
    # calculate Bete_{t} and Alpha_{t} [Y_{it}(0)=beta_{t}X_{it}+alpha_{t}Y_{it}(0)]
    p=length(var_x)
    paramt = data.frame(matrix(0, nrow = T, ncol=p+1))
    colnames(paramt) = c('delt_t',var_x)
    for(t in 1:T0){
      data_temp = data%>%filter(hour_eps==t)%>%dplyr::select(c(var_x, var_y))
      res_ols=lm(pm25~.,data = data_temp)
      paramt[t,]=res_ols$coefficients
    }
    for(t in (T0+1):T){
      data_temp = data_control%>%filter(hour_eps==t)%>%dplyr::select(c(var_x, var_y))
      res_ols=lm(pm25~.,data = data_temp)
      paramt[t,]=res_ols$coefficients
    }
    weight_sum=abs(paramt)[,-1]
  }
  # iterate on t_begin:T, t_begin=1 by default
  T_matched=0
  for(t in t_begin:T){
    # prepare the episode matrix Z1 Z0
    {
      ## X_exo_1: p*1 the covariates of treatment episode at time t
      {
        X_exo_1=data_treated[data_treated$hour_eps==t,var_x_exo]
        rownames(X_exo_1)='treat'
      }
      ## X_en_1: the response for treated's estimated control at lag1 
      {
        if(t<=(T0+1)){
          X_en_1=data.frame(data_treated[data_treated$hour_eps==t, var_x_en])
        }else{
          X_en_1=data.frame(data_treated[data_treated$hour_eps==t, var_x_en])
          colnames(X_en_1)=var_x_en
          X_en_1[1,]=NA
          for(l in 1:K_lag){
            if(t-l<=T0){
              X_en_1[var_x_en[l]]=data_treated[data_treated$hour_eps==t, var_x_en[l]]
            }else{
              X_en_1[var_x_en[l]]=var_y_pre[t-l] #only when the endogenous variable is the lagged response variable
            }
          }
        }
        rownames(X_en_1)='treat'
        colnames(X_en_1)= var_x_en
      }
      ## select the records with complete variables
      {
        data_temp = data_control%>%filter(hour_eps==t)
        id_eps_control_complete = data_temp$id_eps[apply(is.na(data_temp),1,sum)==0]
        row_control = (data_control$hour_eps==t)&(data_control$id_eps%in%id_eps_control_complete)
      }
      ## X_exo_0: p*N the covariates of control episode at time t-1 and t
      {
        X_exo_0=data_control[row_control,var_x_exo]
        rownames(X_exo_0)=data_control[row_control,'id_eps']
      }
      ## X_en_0: the response for control at lag1
      {
        X_en_0=data.frame(data_control[row_control, var_x_en])
        rownames(X_en_0)=1:nrow(X_en_0)
        colnames(X_en_0)=var_x_en
      }
      ## Z1 and Z0
      {
        Z1=t(data.matrix(cbind(X_en_1, X_exo_1)))# (p+1+r_opt)*1
        Z0=t(data.matrix(cbind(X_en_0, X_exo_0)))# (p+1+r_opt)*N
        big.dataframe <- cbind(Z1, Z0)
        std_err_Z <- sqrt(apply(big.dataframe, 1, var))
      }
    }
    # determine the weight
    {
      nvarsV = nrow(Z1)
      if(!is.null(weight_sum)){
        V0=diag(weight_sum[t,], nrow=nvarsV, ncol=nvarsV)
      }else{
        V0=diag(1, nrow=nvarsV, ncol=nvarsV)
      }
      V=V0
      #P=diag(1,nrow=nvarsV,ncol=nvarsV)-rep(1,nvarsV)%*%t(rep(1,nvarsV))/nvarsV
    }
    # solve weighting vector solution.W under the linear constraint
    {
      P=rep(1,ncol(Z0))%*%t(rep(1,ncol(Z0)))
      H <- t(Z0) %*% V %*% (Z0)+P
      c <- -1 * c(t(Z1) %*% V %*% (Z0))
      A <- t(rep(1, length(c)))
      b <- bound
      u <- rep(1, length(c))
      res <- LowRankQP(Vmat = H, dvec = c, Amat = A, bvec = b,uvec = u, method = "LU")
      solution.w <- as.matrix(res$alpha)
    }
    # if the linear constraint can be satisfied, i.e. ,Z1 is within the convex hull of Z0
    {
      Z_diff=Z1-Z0%*%solution.w
      if(mean(abs(Z_diff))<=0.000001){
        solution.w.str=solution.w
        objfun=function(x){
          return(-prod(x))
        }
        Aeq=as.matrix(rbind(Z0,rep(1,dim(Z0)[2])))
        Beq=as.matrix(c(Z1,1))
        lb=rep(0,dim(Z0)[2])
        ub=rep(1,dim(Z0)[2])
        res_nl=NULL
        res_nl=solnl(X=solution.w, objfun=objfun, lb=lb, ub=ub, Aeq=Aeq, Beq=Beq)
        if(!is.null(res_nl)){
          solution.w=res_nl$par
          T_matched=T_matched+1
        }
      }
    }
    # save the results
    {
      if(t==t_begin){
        weight_mat[weight_mat$id_eps_control%in%id_eps_control_complete,1+t]=solution.w
        var_Z_mat=Z0%*%solution.w
        var_diff_mat=Z1-Z0%*%solution.w
        colnames(var_Z_mat)=t
        colnames(var_diff_mat)=t
      }else{
        weight_mat[weight_mat$id_eps_control%in%id_eps_control_complete,1+t]=solution.w
        #colnames(weight_mat)[ncol(weight_mat)]=t
        var_Z_mat=cbind(var_Z_mat, Z0%*%solution.w)
        var_diff_mat=cbind(var_diff_mat, Z1-Z0%*%solution.w)
        colnames(var_diff_mat)[ncol(var_diff_mat)]=t
        colnames(var_Z_mat)[ncol(var_Z_mat)]=t
      }
      var_y_ctl=data_control[row_control, var_y]
      var_y_pre[t]=sum(var_y_ctl*solution.w)
    }
  }  
  # return the results
  return(list(Y0_pre=var_y_pre, weighting_matrix=weight_mat, X0_match=var_Z_mat, X0_match_diff=var_diff_mat,T_matched=T_matched))
}
