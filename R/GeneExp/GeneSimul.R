library(CVST); library(kernrank); library(datamicroarray); library(caret)
source('SVM.R'); source('KLR.R'); source('M_AllKendall.R'); source('GeneSimulate.R');
i=1
iter=1
set.seed(iter)
string = c('alon', 'gordon', 'pomeroy', 'singh', 'petricoin')
datasets = string[i] 
Cancer_result = GeneSimulate(datasets, rep=20, prop=0.7, kernel = c('M_AllKendall', 'kendall_total'), 
                             bandwidth = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1), 
                             C = c(0.01, 0.1, 1, 10, 100, 1000), 
                             lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10), 
                             tol = 10e-6, maxiter = 100, rank = T)

save(Cancer_result, file = paste0('cancer_result_', string[i], '_', iter, '.Rdata' ))
