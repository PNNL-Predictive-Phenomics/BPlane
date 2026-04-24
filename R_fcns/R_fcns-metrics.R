calc_pplus_pmin = function(predictor,response){
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  ones_index = which(response==1)
  zeroes_index = which(response==0)
  
  predictor_ones =  predictor[ones_index]
  predictor_zeroes = predictor[zeroes_index]
  
  p_plus = mean(predictor_ones)
  p_min = mean(predictor_zeroes)
  
  return(list(p_plus=p_plus,p_minus=p_min))
}

# from "metric_functions.R" at 
# https://github.com/lucasvogels33/Review-paper-Bayesian-Structure-Learning-in-GGMs/blob/main/Section%204%20-%20Empirical%20comparison/metric_functions.R
calc_AUC_ROC = function(predictor,response){
  
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #order the vectors so that the predictor is increasing
  predictor.order = order(predictor,decreasing=FALSE)
  predictor.sorted = predictor[predictor.order]
  response.sorted = response[predictor.order]
  
  #determine amount of zeroes and ones
  ones = sum(response)
  zeroes = length(response)-ones
  
  #if there are duplicates
  if (sum(duplicated(predictor.sorted))>0){
    #create a vector with one index for every group of duplicates
    dup_index = cumsum(duplicated(predictor.sorted)==0)
    
    #create a vector sum_vec that sums the true positives in each group of duplicates   
    df <- data.frame(duplicates=dup_index,response.sorted=response.sorted)
    sum_vec = aggregate(response.sorted ~ duplicates, data=df, sum)[,2]
    
    #create a vector that averages the maximum amount of false positives of the current group with the previous group
    fp = cumsum(response.sorted==0)
    df <- data.frame(duplicates=dup_index,fp=fp)
    max_vec = aggregate(fp ~ duplicates, data=df, max)[,2]
    top = c(0,max_vec)
    bottom = c(max_vec,0)
    average_vec = head((top+bottom)/2,-1)
    
    #AUC is the dot product of the two vectors divided by the normalizing constant
    AUC = (sum_vec%*%average_vec)/(ones*zeroes)
    
  }
  
  #if there are no duplicates
  if (sum(duplicated(predictor.sorted))==0){
    fp = cumsum(response.sorted==0)
    AUC = sum(fp * response.sorted)
    AUC = AUC/(zeroes*ones)  
    
  }
  
  return(as.double(AUC))
}

calc_confusion_prob <- function(adj, omega)
{
	conf <- rep(0, 4)
	det <- (omega > .5) |> ifelse(1, 0)
	conf[1] <- (det == 1 & adj == 1 & upper.tri(adj)) |> sum() #tp
	conf[2] <- (det == 1 & adj == 0 & upper.tri(adj)) |> sum() #fp
	conf[3] <- (det == 0 & adj == 0 & upper.tri(adj)) |> sum() #tn
	conf[4] <- (det == 0 & adj == 1 & upper.tri(adj)) |> sum() #fn
	conf
}

calc_metrics <- function(conf)
{
  tp <- conf[1]
  fp <- conf[2]
  tn <- conf[3]
  fn <- conf[4]
  
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  f1 <- 2 * tp / (2*tp + fn + fp)
  mcc <- (tp*tn - fn*fp) / sqrt((tp+fn)*(fp+tn)*(tp+fp)*(fn+tn))
  
  list(tpr=as.double(tpr), fpr=as.double(fpr), sens=as.double(sens), 
       spec=as.double(spec), f1=as.double(f1), mcc=as.double(mcc))
}