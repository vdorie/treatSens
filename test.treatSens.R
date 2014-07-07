data(lalonde)

re74per <- with(lalonde,re74/1000)
re75per <- with(lalonde,re75/1000)
X <- with(lalonde,cbind(married,age,black,hisp,educ,reo74,reo75))#,re75per))
Z <- re75per #with(lalonde,t)
Y <- with(lalonde,re78/1000)

test.norm <- treatSens(Y~Z+X, sensParam = "cor", grid.dim = c(10,10), standardize = T,
                      trt.family = gaussian,
                      resp.family = gaussian,
                      verbose = T, 
                      nsim = 5, data = lalonde)
