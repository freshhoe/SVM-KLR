# SVM (with ksvm) ----
# Cross validation for tuning ----

cv_svm = function(data, params, fold=5, verbose=F) {
  stopifnot(class(data) == "CVST.data" && class(params) == "CVST.params")
  nParams = length(params)
  dimnames = list(as.character(1:fold), names(params))
  results = matrix(0, fold, nParams, dimnames=dimnames)
  
  validationIndex = createFolds(data$y, k = fold)
  for (ind in 1:nParams) {
    p = params[[ind]]
    if (p$kernel == "M_AllKendall"){
      kernel = eval(parse(text = p$kernel))
      curK <- as.kernelMatrix(kernel(bandwidth = p$bandwidth, data$x))
    } else if (p$kernel == "kendall_total") {
      kernel = eval(parse(text = p$kernel))
      curK <- as.kernelMatrix(kernel(data$x))
    } else {
      kernel = eval(parse(text = p$kernel))
      curK <- as.kernelMatrix(kernel(data$x))
    }
    for (f in 1:fold) {
      curTrain = getSubset(data, -validationIndex[[f]])
      curTest = getSubset(data, validationIndex[[f]])
      curTrainK = as.kernelMatrix(curK[-validationIndex[[f]],-validationIndex[[f]]])
      curmodel = ksvm(curTrainK, data$y[-validationIndex[[f]]], kernel='matrix', type = "C-svc", C = p$C, scaled = F)
      if (p$kernel == "M_AllKendall") {
        curTestK = as.kernelMatrix(M_AllKendall(bandwidth = p$bandwidth, curTest$x, curTrain$x))
        curpred  = predict(curmodel, curTestK)
      } else if (p$kernel == "kendall_total") {
        curTestK = as.kernelMatrix(curK[validationIndex[[f]], -validationIndex[[f]]][,SVindex(curmodel), drop = F])
        curpred  = predict(curmodel, curTestK)
      }
      results[f, ind] = mean(curTest$y != curpred)
    }
    if (verbose) {
      cat(names(params)[ind], "(", mean(results[, ind]), ")\n")
    }
  }
  winner = which.min(apply(results, 2, mean))
  if (length(winner) == 0) {
    return(NULL)
  } else {
    return(params[winner])
  }
}