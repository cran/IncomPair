library(testthat)
library(IncomPair)

n=20
n1=15
n2=10
r=0.8
xp=rnorm(n)
yp=r*xp+(1-r)*rnorm(n)
xu=rnorm(n1)
yu=rnorm(n2)
mu=0

a1<-permb(xp,yp,xu,yu,r,mu,method="EH",alternative="two.sided",verbose = TRUE)
a2<-parmb(xp,yp,xu,yu,r,mu,method="Zb",alternative="two.sided",verbose = TRUE)
a3<-rankb(xp,yp,xu,yu,mu,method="Rankw",alternative="two.sided",verbose = TRUE)

test_that("p-values within range",{
  expect_equal(a1@Pval,.5,tolerance=.5)
  expect_equal(a2@Pval,.5,tolerance=.5)
  expect_equal(a3@Pval,.5,tolerance=.5)
})
