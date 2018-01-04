library(CVST); library(kernrank); library(MASS); library(caret)
source('SVM.R'); source('KLR.R'); source('M_AllKendall.R'); source('simul_func.R');
i=1
j=1
k=1
iter=1
set.seed(iter)
n = c(100, 10000); p = c(10, 50, 250); D = c(1, 1.5, 2, 2.5, 2.75, 3); rho = c(0, 0.5, -0.5); r = c(0.1, 0.3, 0.5); 
simul_dat = list(train = simul_data(n[1], p[i], rho[j], D[6], r[k]), test = simul_data(n[2], p[i], rho[j], D[6], r[k])) 

simul_result = simul_func(simul_dat, rep = 10, kernel = c('M_AllKendall', 'kendall_total'), 
                          bandwidth = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1), 
                          C = c(0.01, 0.1, 1, 10, 100, 1000, 10000), 
                          lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), 
                          tol = 10e-6, maxiter = 100, rank = T)

save(simul_result, file = paste0('result_', p[i], '_', rho[j], '_', r[k], '_', iter, '.Rdata' ))