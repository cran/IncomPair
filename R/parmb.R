#' @title Parametric Based Tests for Incomplete Paired Data
#'
#' @description The function performs testing the Hypothesis of equality of means for the incomplete pairs setting data. The function uses six test statistics that were proposed for testing the equality of the means of a bivariate normal distribution with unknown common variance and correlation coefficient when observations are missing on both variates. These function includes Lin and Stivers (1974, Ts), Bhoj (1989, pp. 282, Z), Bhoj (1989, pp. 282, Zb),  Bhoj (1989, pp. 283, T), Bhoj (1989, pp. 283, Zh) and  Bhoj (1989, pp 284, Zls). For more details, information of the functions see Bhoj (1989).
#'
#' @param xp,yp (non-empty) numeric vectors of data values of the the complete pairs
#' @param xu a numeric vector of data on x only
#' @param yu a numeric vector of data on y only
#' @param r a number indicating the correlation between the complete pairs
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test
#' @param method a character string specifying the different type of methods, must be one of "Zb" (default), "Zb","T","Tls" ,"Zls","Zh"
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' @param verbose if TRUE, show the test used, null and alternative hypotheses in addition to the p-value
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
#' parmb(xp,yp,xu,yu,r,mu,method="Zb",alternative="two.sided", verbose=TRUE)
#'
#' @references
#' Bhoj, D. S. (1991). Testing equality of means in the presence of correlation and missing data. Biometrical journal, 33(1), 63-72.
#' Ekbohm, G. (1976). On comparing means in the paired case with incomplete data on both responses. Biometrika 63(2), 299-304.
#' Lin, P. E., & Stivers, L. E. (1974). On difference of means with incomplete data. Biometrika, 61(2), 325-334.
#'
#' @export parmb

parmb <- function(xp, yp, xu, yu, r, mu=NULL, method="Zb",
                  alternative="two.sided", verbose=TRUE)


{

  met <- pmatch(method, c("Z", "Zb", "T", "Tls", "Zls", "Zh" ))
  if (is.na(met)) stop("\nmethod MUST be one of:'Z', 'Zb', 'T', 'Tls', 'Zls', or 'Zh' \n")

  if (missing(r)) {
    r <- cor(xp,yp)
  }

  else {
    if ((length(xp) <= 1) |(length(yp) <= 1)) {
      stop("\n There are no paired observations\n")
    }
  }

  alt <- pmatch(alternative, c("two.sided","less", "greater"))
  if (is.na(alt)) stop("\nalternative MUST be one of:'two.sided', 'less', or 'greater' \n")



  n=length(xp)
  n1=length(xu)
  n2=length(yu)

  #Put data together
  x=c(xp,xu)
  y=c(yp,yu)


  #Components for Statistics
  nx=sum(xp)
  ny=sum(yp)
  n1x=sum(xu)
  n2y=sum(yu)
  nn1x=sum(x)
  nn2y=sum(y)
  a11=sum((xp-mean(xp))^2)
  a22=sum((yp-mean(yp))^2)
  a12=sum((xp-mean(xp))*(yp-mean(yp)))
  b1=sum((xu-mean(xu))^2)
  b2=sum((yu-mean(yu))^2)
  c1=sum((x-mean(x))^2)
  c2=sum((y-mean(y))^2)
  r=a12/(sqrt(a11*a22))
  u=(2*a12)/(a11+a22)

  t1num = (mean(xp)-mean(yp))*(sqrt(n))
  t1denom = ((a11+a22-(2*a12))/(n-1))^(1/2)
  t1=t1num/t1denom
  f1= n-1

  t2num= (mean(xu)-mean(yu))
  t2denom= (((b1+b2)/((n1+n2-2)))*((1/n1)+(1/n2)))^(.5)
  t2=t2num/t2denom
  f2= n1 + n2 -2

  dhat2= 2*mean(xu) - mean(xp) - mean(yp)
  dhat3= mean(xp) + mean(yp) - 2*mean(yu)
  w= n1*(n+(1+ r)*n2)*((n*(n1+n2)+(2*(1+r)*n1*n2))^-1)
  dhat = (w*dhat2)+((1-w)*dhat3)
  s= (1+u)/2
  t3denom1=((4*s*(b1+b2)+a11+a22+2*a12)/(n+n1+n2-3))
  t3denom2=((w^2)/(s*n1))+(((1-w)^2)/(s*n2))+(((1-2*w)^2)/(n))
  t3=dhat/((t3denom1*t3denom2)^(1/2))
  f3= n+n1+n2-3


  m1=n/(n+n1)
  m2=n/(n+n2)



  #Test Statistics for paper A

  if(met==1)

  {
    #Z
    lambda = (1 + sqrt((2*n1*n2*(1-r))/(n*(n1+n2))))^-1
    F1= 1 + ((2*t1^2)/f1) + (2*t1/sqrt(f1))*sqrt(1+((t1^2)/f1))
    F2= 1 + ((2*t2^2)/f2) + (2*t2/sqrt(f2))*sqrt(1+((t2^2)/f2))
    F3= 1 + ((2*t3^2)/f3) + (2*t3/sqrt(f3))*sqrt(1+((t3^2)/f3))
    U1= (1-(2/(9*f1)))*((F1^(1/3))-1)*((2/(9*f1))*((F1^(2/3))+1))^(-1/2)
    U2= (1-(2/(9*f2)))*((F2^(1/3))-1)*((2/(9*f2))*((F2^(2/3)+1)))^(-1/2)
    U3= (1-(2/(9*f3)))*((F3^(1/3))-1)*((2/(9*f3))*((F3^(2/3)+1)))^(-1/2)
    Z= (lambda*U1+(1-lambda)*U2)*(lambda^2+(1-lambda)^2)^(-1/2)


    if(alt==1)  {pvalue=2*(1-pnorm(abs(Z)))}

    if (alt==2) {pvalue=pnorm(Z)}

    if (alt==3) {pvalue=1-pnorm(Z)}
  }


  if(met==2)
  {
    #Zb
    F1= 1 + ((2*t1^2)/f1) + (2*t1/sqrt(f1))*sqrt(1+((t1^2)/f1))
    F2= 1 + ((2*t2^2)/f2) + (2*t2/sqrt(f2))*sqrt(1+((t2^2)/f2))
    F3= 1 + ((2*t3^2)/f3) + (2*t3/sqrt(f3))*sqrt(1+((t3^2)/f3))
    U1= (1-(2/(9*f1)))*((F1^(1/3))-1)*((2/(9*f1))*((F1^(2/3))+1))^(-1/2)
    U2= (1-(2/(9*f2)))*((F2^(1/3))-1)*((2/(9*f2))*((F2^(2/3)+1)))^(-1/2)
    U3= (1-(2/(9*f3)))*((F3^(1/3))-1)*((2/(9*f3))*((F3^(2/3)+1)))^(-1/2)
    lambdab= (1+sqrt((n1*n2*(1-r))/((2*n*n2*w^2)+(2*n*n1*(1-w)^2)+(n1*n2*(1-2*w)^2)*(1+r))))^-1
    Zb= (lambdab*U1+(1-lambdab)*U3)*(lambdab^2+(1-lambdab)^2)^(-1/2)


    if(alt==1)  {pvalue=2*(1-pnorm(abs(Zb)))}

    if (alt==2) {pvalue=pnorm(Zb)}

    if (alt==3) {pvalue=1-pnorm(Zb)}

  }


  if(met==3)
  {
    #T
    m1=n/(n+n1)
    m2=n/(n+n2)
    xy= m1*(mean(xp)-mean(xu))+mean(xu)-(m2*(mean(yp)-mean(yu))+mean(yu))
    k1=((m1^2)+(m2^2)-(2*m1*m2*u))/((1-m1)^2)
    k2=((m1^2)+(m2^2)-(2*m1*m2*u))/((1-m2)^2)
    Tdenom1=(((m1^2)*a11+(m2^2)*a22-(2*m1*m2*a12))+(k1*((1-m1)^2)*b1)+(k2*((1-m2)^2)*b2))/(n+n1+n2-3)
    Tdenom2=((1/n)+(1/(k1*n1))+(1/(k2*n2)))
    T=(xy)/((Tdenom1*Tdenom2)^(1/2))
    f4=n+n1+n2-3


    if(alt==1)  {pvalue=2*(1-pt(abs(T),f4))}

    if (alt==2) {pvalue=pt(T,f4)}

    if (alt==3) {pvalue=1-pt(T,f4)}

  }


  if(met==4)
  {

    #Tls
    m1=n/(n+n1)
    m2=n/(n+n2)
    xy= m1*(mean(xp)-mean(xu))+mean(xu)-(m2*(mean(yp)-mean(yu))+mean(yu))
    Tlsdenom1=((1/(n+n1))+(1/(n+n2))-((2*n*r)/((n+n1)*(n+n2))))
    Tlsdenom2=((c1+b2)/(n+n1+n2-2))
    Tls=(xy)/((Tlsdenom1*Tlsdenom2)^(1/2))
    f4=n+n1+n2-4


    if(alt==1)  {pvalue=2*(1-pt(abs(Tls),f4))}

    if (alt==2) {pvalue=pt(Tls,f4)}

    if (alt==3) {pvalue=1-pt(Tls,f4)}

  }


  if(met==5)
  {

    #Zh
    a=(m1+u*f2*(1-f1))*(1-(u^2)*(1-m1)*(1-m2))^-1
    b=(m2+u*f1*(1-f2))*(1-(u^2)*(1-m1)*(1-m2))^-1
    s2hat=(a11+a22+b1+b2)*(2*n+n1+n2-4)^-1
    zhnum=(a*mean(xp))+(1-a)*mean(xu)-b*mean(yp)-(1-b)*mean(yu)
    zhdenom=(s2hat*((a^2)+m1*(((1-a)^2)/(1-m1))+b^2+m2*(((1-b)^2)/(1-m2))-2*a*b*u)/n)^(1/2)
    Zh=zhnum/zhdenom


    if(alt==1)  {pvalue=2*(1-pt(abs(Zh),(n-1)))}

    if (alt==2) {pvalue=pt(Zh,(n-1))}

    if (alt==3) {pvalue=1-pt(Zh,(n-1))}

  }


  if(met==6)

  {

    #Zls
    g=n*(n+n2+(n1*a12)/a11)/((n+n1)*(n+n2)-n1*n2*r^2)
    h=n*(n+n1+(n2*a12)/a22)/((n+n1)*(n+n2)-n1*n2*r^2)
    V1=((((g^2)/n)+(((1-g)^2)/n1))*a11+(((h^2)/n)+(((1-h)^2)/n2))*a22-(2*h*g*a12)/n)/f1
    Zls=(g*mean(xp)+(1-g)*mean(xu)-h*mean(yp)-(1-h)*mean(yu))/(sqrt(V1))


    if(alt==1)  {pvalue=2*(1-pt(abs(Zls),n))}

    if (alt==2) {pvalue=pt(Zls,n)}

    if (alt==3) {pvalue=1-pt(Zls,n)}

  }

  setClass("ans",
           representation(
             Title = "character",
             Nhypothesis = "character",
             Ahypothesis = "character",
             Pval = "numeric"),
           where = topenv(parent.frame())
  )

  if (method == "Z"){
    display<-new("ans",
                 Title = "Bhoj (1989, pp. 282), Equ. 2.4, Z",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))

  }

  if (method == "Zb"){
    display<-new("ans",
                 Title = "Bhoj (1989, pp. 282), Equ. 2.7, Zb",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))
  }

  if (method == "T"){
    display<-new("ans",
                 Title = "Bhoj (1989, pp. 283), T",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))
  }

  if (method == "Tls"){
    display<-new("ans",
                 Title = "Bhoj (1989, pp. 283), Tls",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))
  }

  if (method == "Zh"){
    display<-new("ans",
                 Title = "Bhoj (1989, pp. 283), Zh",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))
  }

  if (method == "Zls"){
    display<-new("ans",
                 Title = "Bhoj (1989, pp. 284), Zls",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))
  }

  if(verbose==TRUE) {
    setMethod("show", "ans",
            function(object){
              cat("\n", object@Title, "\n", "\n",
                  object@Nhypothesis, "\n",
                  object@Ahypothesis, "\n", "\n",
                  "p-value =", object@Pval, "\n","\n")
            }, where = topenv(parent.frame()))
    display} else {cat("p-value =",display@Pval,"\n")}


}
