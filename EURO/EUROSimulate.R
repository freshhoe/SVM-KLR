# EURO barometer data simulation function ----

EUROSimulate = function (data, rep, prop, kernel, bandwidth, C, lambda, tol, maxiter, rank = F){
  Start = paste("Starting the", rep, "times iteration :", date(), "\n")
  cat(Start)
  call <- match.call()
  
  # Data transformation ----
  if(rank == T) data$train$x = floor(t(apply(data$x, 1, rank)))
  
  # Mallows result ----
  M_piResult = list(); M_pred_klr = list(); M_pred_svm = list(); 
  M_MISMat_klr = matrix(0, rep, 1); M_MISMat_svm = matrix(0, rep, 1)
  
  # Kendall result ----
  K_piResult = list(); K_pred_klr = list(); K_pred_svm = list(); 
  K_MISMat_klr = matrix(0, rep, 1); K_MISMat_svm = matrix(0, rep, 1)
  
  # CVST.params for Mallows ----
  M_klr_par = constructParams(kernel = kernel[1], bandwidth = bandwidth, lambda=lambda, tol=tol, maxiter=maxiter) 
  M_svm_par = constructParams(kernel = kernel[1], bandwidth = bandwidth, C = C)  
  
  # CVST.params for Kendall ----
  K_klr_par = constructParams(kernel = kernel[2], lambda=lambda, tol=tol, maxiter=maxiter) 
  K_svm_par = constructParams(kernel = kernel[2], C = C) 
  
  for (i in 1:rep){
    train_set = createDataPartition(data$y, p = prop, list = T)
    trainX = data$x[train_set[[1]],]
    testX = data$x[-train_set[[1]],]
    trainY = data$y[train_set[[1]]]
    testY = data$y[-train_set[[1]]]
    Train = constructData(trainX, trainY)
    Test = constructData(testX, testY)
    
    #.. Mallows tuning ----
    M_KLR = paste("Mallows KLR tuning...:", date(), "\n"); cat(M_KLR)
    M_opt_klr = cv_klr(Train, M_klr_par) 
    M_SVM = paste("Mallows SVM tuning...:", date(), "\n"); cat(M_SVM)
    M_opt_svm = cv_svm(Train, M_svm_par)
    
    M_kernel_func = eval(parse(text = kernel[1]))
    M_klrK <- as.kernelMatrix(M_kernel_func(bandwidth = M_opt_klr[[1]]$bandwidth, trainX))
    M_svmK <- as.kernelMatrix(M_kernel_func(bandwidth = M_opt_svm[[1]]$bandwidth, trainX))
    
    #..... KLR Prediction ----
    M_klr_TrainK = Matrix(M_klrK)
    M_klr_model = Klogreg(trainX, kernel[1], M_klr_TrainK, trainY, 
                          M_opt_klr[[1]]$lambda, M_opt_klr[[1]]$tol, M_opt_klr[[1]]$maxiter)
    M_klr_TestK = M_kernel_func(bandwidth = M_opt_klr[[1]]$bandwidth, testX, trainX)
    M_pred_klr[[i]] = Klogreg_predict(M_klr_model, M_klr_TestK)
    M_piResult[[i]] = M_pred_klr[[i]]$pi
    M_MISMat_klr[i,] = sum(M_pred_klr[[i]]$result != testY) / length(testY)
    
    #..... SVM Prediction ----
    M_svm_TrainK = M_svmK
    M_svm_model = ksvm(M_svm_TrainK, trainY, kernel='matrix', type = "C-svc", C = M_opt_svm[[1]]$C, scaled = F)
    M_svm_TestK = as.kernelMatrix(M_kernel_func(bandwidth = M_opt_svm[[1]]$bandwidth, testX, trainX))
    M_pred_svm[[i]] = predict(M_svm_model, M_svm_TestK)
    M_MISMat_svm[i,] = sum(M_pred_svm[[i]] != testY) / length(testY)
    
    #.. Kendall tuning ----
    K_KLR = paste("Kendall KLR tuning...:", date(), "\n"); cat(K_KLR)
    K_opt_klr = cv_klr(Train, K_klr_par) 
    K_SVM = paste("Kendall SVM tuning...:", date(), "\n"); cat(K_SVM)
    K_opt_svm = cv_svm(Train, K_svm_par) 
    
    K_kernel_func = eval(parse(text = kernel[2]))
    K_klrK <- as.kernelMatrix(K_kernel_func(trainX))
    K_svmK <- K_klrK
    
    #..... KLR Prediction ----
    K_klr_TrainK = Matrix(K_klrK)
    K_klr_model = Klogreg(trainX, kernel[2], K_klr_TrainK, trainY, 
                          K_opt_klr[[1]]$lambda, K_opt_klr[[1]]$tol, K_opt_klr[[1]]$maxiter)
    K_klr_TestK = Test
    K_pred_klr[[i]] = Klogreg_predict(K_klr_model, Test)
    K_piResult[[i]] = K_pred_klr[[i]]$pi
    K_MISMat_klr[i,] = sum(K_pred_klr[[i]]$result != testY) / length(testY)
    
    #..... SVM Prediction ----
    dat = data.frame(y = trainY, trainX) 
    kernel_func = eval(parse(text = kernel[2]))
    class(kernel_func) = 'kernel'
    K_svm_model = ksvm(dat$y ~., data = dat, kernel = kernel_func, type = "C-svc", C = K_opt_svm[[1]]$C, scaled = F)
    K_pred_svm[[i]]  = predict(K_svm_model, testX)
    K_MISMat_svm[i,] = sum(K_pred_svm[[i]] != testY) / length(testY)
    
    ith.message = paste(i, "result calculation is complete\n")
    cat(ith.message)
  }
  
  Close = paste("Finishing the", rep, "times iteration :", date(), "\n")
  cat(Close)
  StartClose = c(Start, Close)
  
  # Mallows result ----
  # .. KLR Result ----
  M_Average_klr = round(colMeans(M_MISMat_klr), 4); M_SE_klr = round(sd(M_MISMat_klr[,1])/sqrt(nrow(M_MISMat_klr)), 4)
  M_MeanSE_klr  = paste0(M_Average_klr, ' (', M_SE_klr, ')')
  M_Output_klr  = list(M_pred_klr = M_pred_klr, M_MISMat_klr = M_MISMat_klr, M_piResult = M_piResult, M_MeanSE_klr = M_MeanSE_klr)
  
  # .. SVM Result ----
  M_Average_svm = round(colMeans(M_MISMat_svm), 4); M_SE_svm = round(sd(M_MISMat_svm[,1])/sqrt(nrow(M_MISMat_svm)), 4)
  M_MeanSE_svm  = paste0(M_Average_svm, ' (', M_SE_svm, ')')
  M_Output_svm  = list(M_pred_svm = M_pred_svm, M_MISMat_svm = M_MISMat_svm, M_MeanSE_svm = M_MeanSE_svm)
  
  # Kendall result ----
  # .. KLR Result ----
  K_Average_klr = round(colMeans(K_MISMat_klr), 4); K_SE_klr = round(sd(K_MISMat_klr[,1])/sqrt(nrow(K_MISMat_klr)), 4)
  K_MeanSE_klr  = paste0(K_Average_klr, ' (', K_SE_klr, ')')
  K_Output_klr  = list(K_pred_klr = K_pred_klr, K_MISMat_klr = K_MISMat_klr, K_piResult = K_piResult, K_MeanSE_klr = K_MeanSE_klr)
  
  # .. SVM Result ----
  K_Average_svm = round(colMeans(K_MISMat_svm), 4); K_SE_svm = round(sd(K_MISMat_svm[,1])/sqrt(nrow(K_MISMat_svm)), 4)
  K_MeanSE_svm  = paste0(K_Average_svm, ' (', K_SE_svm, ')')
  K_Output_svm  = list(K_pred_svm = K_pred_svm, K_MISMat_svm = K_MISMat_svm, K_MeanSE_svm = K_MeanSE_svm)
  
  Total_Output = list(StartClose = StartClose, call = call, K_Output_klr = K_Output_klr, K_Output_svm = K_Output_svm,
                      M_Output_klr = M_Output_klr, M_Output_svm = M_Output_svm)
  return(Total_Output)
}