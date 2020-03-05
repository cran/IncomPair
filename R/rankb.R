#' @title Rank Based Tests for Incomplete Paired Data
#'
#' @description The function performs testing the hypothesis of equality of means for the incomplete pairs setting data. The function uses a rank-based procedure for parameter estimation and hypothesis testing when the data are a mixture of paired observations and independent samples. The rank-based methods combine Wilcoxon signed-rank statistics and Wilcoxon-Mann-Whitney two-sample procedures. These methods were developed by Dubnicka, Blair and Hettmansperger (2002).
#'
#' @param xp,yp (non-empty) numeric vectors of data values of the the complete pairs
#' @param xu a numeric vector of data on x only
#' @param yu a numeric vector of data on y only
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test
#' @param method a character string specifying the different type of methods, must be one of "Ranku" or "Rankw"
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @param verbose if TRUE, show the test used, null and alternative hypotheses in addition to the p-value
#'
#'
#' @return NULL
#'
#' @examples
#' n=20
#' n1=15
#' n2=10
#' r=0.8
#' xp=rnorm(n)
#' yp=r*xp+(1-r)*(rnorm(n))
#' xu=rnorm(n1)
#' yu=rnorm(n2)
#' mu=0
#' rankb(xp,yp,xu,yu,mu,method="Rankw",alternative="two.sided",verbose=TRUE)
#'
#' @references
#' Dubnicka, S. R., Blair, R. C., & Hettmansperger, T. P. (2002). Rank-based procedures for mixed paired and two-sample designs. Journal of Modern Applied Statistical Methods, 1(1), 6.
#'
#' @export rankb

rankb <- function(xp, yp, xu, yu, mu=NULL, method="Ranku",
                  alternative="two.sided", verbose=TRUE)

{
  met <- pmatch(method, c("Ranku","Rankw"))
  if (is.na(met))
    stop("\nmethod MUST be one of: 'Ranku', or 'Rankw'\n")
  alt <- pmatch(alternative, c("two.sided","less", "greater"))
  if (is.na(alt))
    stop("\nalternative MUST be one of: 'two.sided', 'less', or 'greater'\n")

  if ((length(xp) <= 1) |(length(yp) <= 1)) {
    stop("\nThere are no paired observations\n")
  }



  n1=length(xu)
  n2=length(yu)

  #### Rank test for both weighted and unweighetd tests

  D=xp-yp
  n=length(D)
  R0=rank(abs(D))
  xy=c(xu,yu)
  Ru=rank(xy)

  R1=sum(Ru[1:n1])

  Tp0=0
  for (ik in 1:n)
  {

    if(D[ik]> 0) Tp0=Tp0+R0[ik]

  }

  Tu=R1-n1*(n1+1)/2

  Tr=Tp0+Tu


  mu1=n*(n+1)/4 +n1*n2/2

  vu=n*(n+1)*(2*n+1)/24 + n1*n2*(n1+n2+1)/12

  TRu=(Tr-mu1)/sqrt(vu)

  N=n1+n2

  aw=2*N/((n*N+2*n1*n2)*(n+1))
  bw=2/(n*N+2*n1*n2)

  a1w=n*N/(n*N+2*n1*n2)
  b1w=2*n1*n2/(n*N+2*n1*n2)

  vw=(aw^2)*(n*(n+1)*(2*n+1))/24 + (bw^2)*(n1*n2*(N+1))/12

  Tw=aw*Tp0 +bw*Tu

  Twe=aw*(n*(n+1))/4 + bw*(n1*n2)/2

  Tw1=(Tw-0.5)/(sqrt(vw))

  if(met==1)
  {

    if(alt==1)  {pvalue=2*(1-pnorm(abs(TRu)))}

    if (alt==2) {pvalue=pnorm(TRu)}

    if (alt==3) {pvalue=1-pnorm(TRu)}

  }


  if(met==2)
  {

    if(alt==1)  {pvalue=2*(1-pnorm(abs(Tw1)))}

    if (alt==2) {pvalue=pnorm(Tw1)}

    if (alt==3) {pvalue=1-pnorm(Tw1)}

  }

  setClass("ans",
           representation(
             Title = "character",
             Nhypothesis = "character",
             Ahypothesis = "character",
             Pval = "numeric"),
           where = topenv(parent.frame())
  )

  if (method == "Ranku"){
    display<-new("ans",
                 Title = "Rank Based Unweighted Test",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))

  }

  if (method == "Rankw"){
    display<-new("ans",
                 Title = "Rank Based Weighted Test",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))
  }

  if(verbose==TRUE){
    setMethod("show", "ans",
            function(object){
              cat("\n", object@Title, "\n", "\n",
                  object@Nhypothesis, "\n",
                  object@Ahypothesis, "\n", "\n",
                  "p-value =", object@Pval, "\n","\n")
            }, where = topenv(parent.frame()))
    display} else {cat("p-value =",display@Pval,"\n")}

}
