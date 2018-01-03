# SVM-KLR
R을 이용한 순위 데이터에 대한 SVM과 KLR의 성능 비교 

## Note
본 결과는 학위논문의 일환이며 `CVST` 패키지와 `kerrank` 패키지를 기반으로 작성된 코드입니다.

## Installation
`kerrank` 패키지 설치
```r
# install.packages("devtools")
devtools::install_github("YunlongJiao/kernrank")
```

## Data
#### Gene expression data
[kendallkernel_demo](https://github.com/YunlongJiao/kendallkernel_demo/tree/master/geneexpr/data) 참조

#### Euro barometer data


## Reference
- Genman, D. et. al (2004). Classifying gene expression profiles from pairwise mRNA comparisons.
- Jiao, Y. and Vert, J.-P. (2015). The Kendall and Mallows kernels for permutations.
- Mania, H. et. al (2016). On kernel methods for covariates that are rankings. 
- Kimeldorf, G. and Wahba, G. (1971). Some Results on Tchebycheffian Spline Functions.
- Zhu, J. and Hastie, T. (2005). Kernel Logistic Regression and Import Vector Machine.
