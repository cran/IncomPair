#' @title Permutation Based Tests for Incomplete Paired Data
#'
#' @description The function Performs testing the hypothesis of equality of means for the incomplete pairs setting data. The function combines the observed mean difference for the complete pairs with the difference between the two means of the independent samples. The function implements two different nonparametric tests based on permutation tests that were proposed by Einsporn and Habtzghi (2013), and Maritz (1995). The two methods are denoted by EH and Maritz, respectively.
#'
#' @param xp,yp (non-empty) numeric vectors of data values of the the complete pairs
#' @param xu a numeric vector of data on x only
#' @param yu a numeric vector of data on y only
#' @param r a number indicating the correlation between the complete pairs
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test
#' @param method a character string specifying the different type of methods, must be one of "EH" (default) or "Maritz"
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
#' yp=r*xp+(1-r)*rnorm(n)
#' xu=rnorm(n1)
#' yu=rnorm(n2)
#' mu=0
#' permb(xp,yp,xu,yu,r,mu,method="EH",alternative="two.sided",verbose=TRUE)
#'
#' @references
#' 1 Einsporn, R. L., & Habtzghi, D. (2013). Combining paired and two-sample data using a permutation test. Journal of Data Science, 11(4), 767-779.
#' 2 Maritz, J. S. (1995). A permutation paired test allowing for missing values. Australian Journal of Statistics, 37(2), 153-159.
#' 3 Woolson, R., Leeper, J., Cole, J. and Clarke, W. (1976). A Monte Carlo investigation of a statistic for a bivariate missing data problem. Communications in Statistics - Theory and Methods A5, 681-688.
#'
#' @export permb

permb <- function(xp, yp, xu, yu, r, mu=NULL, method="EH", alternative="two.sided", verbose=TRUE)

{
  met <- pmatch(method, c("EH", "Maritz"))
  if (is.na(met))
    stop("\nmethod MUST be one of: 'EH', 'Maritz'\n")
  alt <- pmatch(alternative, c("two.sided","less", "greater"))
  if (is.na(alt))
    stop("\nalternative MUST be one of: 'two.sided', 'less', or 'greater' \n")

  if (missing(r)) {
    r <- cor(xp,yp)
  }

  else {
    if ((length(xp) < 1) |(length(yp) < 1)) {
      stop("\n There are no paired observations\n")
    }
  }


  if (met==1)
  {
    # this function computes the mean of the differences for the unpaired sample

    # This part compute permutation test for paired data
    paired <- function(x1,x2)
    {
      n=length(x1)
      y1=numeric(n)
      y2=numeric(n)

      for(i in 1:n)
      { u=runif(1)
      if ((u-0.5)<=0) {
        y1[i]=x1[i]
        y2[i]=x2[i]

      } else {

        if((u-0.5)>0) {

          y1[i]=x2[i]
          y2[i]=x1[i]
        }
      }
      }

      axp=mean(y1-y2)


      return(axp)

    }


    unpaired <- function(x1,x2)

    {

      m1=length(x1)
      m2=length(x2)


      tc=c(x1,x2)
      n=length(tc)


      y=sample(tc,replace=F)
      ax=mean(y[1:m1])
      ay=mean(y[(m1+1):n])

      # the difference between the sample means
      cmean=ax-ay

      return(cmean)

    }

    # weight for equal variance

    permt <- function(xp,yp,xu,yu,r)

    {
      n=length(xp)
      n1=length(xu)
      n2=length(yu)
      w=(1/n2+1/n1)/((2-2*r)/n+(1/n2+1/n1))
      if(n==0) w=0
      if((n1==0)&&(n2==0)) w=1


      Tp=paired(xp,yp)

      Tu=unpaired(xu,yu)

      T0=w*(Tp) +(1-w)*(Tu)

      return(T0)
    }

    np=5000
    T0=replicate(np,permt(xp,yp,xu,yu,r))

    n=length(xp)
    n1=length(xu)
    n2=length(yu)
    w=(1/n2+1/n1)/((2-2*r)/n +(1/n2+1/n1))
    if(n==0) w=0

    if((n1==0)&&(n2==0)) w=1


    # the value of the observed test statistic

    Tobs=w*(mean(xp)-mean(yp))+(1-w)*(mean(xu)-mean(yu))

    result=c(T0,Tobs)



    #The p-value is the chance of seeing a test statistic at least as
    #large as what was observed (.667). This can be estimated by the
    #proportion of permutation statistics that were as
    #large or larger than the observed one

    if(alt==1)  { pvalue=sum(abs(T0)>=abs(Tobs))/length(result)

    }

    if(alt==2)  { pvalue=sum(T0<=Tobs)/length(result)

    }

    if(alt==3)  { pvalue=sum(T0>=Tobs)/length(result)

    }
  }


  if(met==2)
  {


    #  Maritz's Permutation test


    permm <- function(xp,yp,xu,yu)

    {
      n=length(xp)
      n1=length(xu)
      n2=length(yu)


      #
      x=c(xu, rep(0,n2),xp)
      y=c(rep(0,n1),yu,yp)




      m=n1+n2+n

      y1=numeric(m)
      y2=numeric(m)

      for(i in 1:m)

      {
        T2=sample(c(x[i],y[i]),replace=F)
        y1[i]=T2[1]
        y2[i]=T2[2]

      }


      b=length(y1[y1!=0])
      b2=length(y2[y2!=0])

      xbar=sum(y1)/b
      ybar=sum(y2)/b2



      diff=xbar-ybar


      return(diff)

    }




    #Maritz

    np=5000
    tm0=replicate(np,permm(xp,yp,xu,yu))


    # this part computes the observed value of Maritz

    n=length(xp)
    n1=length(xu)
    n2=length(yu)

    xl=c(xu, rep(0,n2),xp)
    yl=c(rep(0,n1),yu,yp)
    bl=c(rep(1,n1),rep(0,n2),rep(1,n))


    #ayp <- numeric(np)

    xbar1=sum(xl*bl)/sum(bl)
    ybar1=(sum(xl)+sum(yl)-sum(xl*bl))/(n1+n2+2*n-sum(bl))


    Dm0=xbar1-ybar1



    resultm=c(tm0,Dm0)


    if(alt==1)  { pvalue=sum(abs(tm0)>=abs(Dm0))/length(resultm)

    }

    if(alt==2)  {pvalue=sum(tm0<=Dm0)/length(resultm)

    }

    if(alt==3)  { pvalue=sum(tm0>=Dm0)/length(resultm)

    }
  }

  setClass("ans",
           representation(
             Title = "character",
             Nhypothesis = "character",
             Ahypothesis = "character",
             Pval = "numeric"),
           where = topenv(parent.frame())
  )

  if (method == "EH"){
    display<-new("ans",
                 Title = "Permutation Based Test Using Einsporn and Habtzghi's Approach",
                 Nhypothesis = "Ho: The true difference between means is 0.",
                 Ahypothesis = switch(alternative, two.sided = "Ha: The true difference between means is not equal to 0.",
                                      less = "Ha: The true difference between means is less than 0.",
                                      greater = "Ha: The true difference between means is greater than 0."),
                 Pval = round(pvalue,digits = 5))

  }

  if (method == "Maritz"){
  display<-new("ans",
               Title = "Permutation Based Test Using Maritz's Approach",
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
    display} else {
      cat("p-value =",display@Pval,"\n")}

}
