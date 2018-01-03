library(CVST); library(kernrank); library(MASS); library(caret)
source('SVM.R'); source('KLR.R'); source('M_AllKendall.R'); source('EUROSimulate.R');
iter=1
set.seed(iter)
load('survey.RData')
EURO_result = EUROSimulate(survey, rep=20, prop=0.02395, kernel = c('M_AllKendall', 'kendall_total'), 
                           bandwidth = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1), 
                           C = c(0.01, 0.1, 1, 10, 100, 1000, 10000), 
                           lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), 
                           tol = 10e-6, maxiter = 100, rank = F)

save(EURO_result, file = paste0('euro_result_', iter, '.Rdata' ))
