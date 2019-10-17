context("Generic tests")
test_that("Sanity checks", {

  var <- array(rnorm(16), c(2, 2, 2, 2))
  names(dim(var)) <- c("memb", "lon", "sdates", "lat")
  clim <- apply(var, c(1, 2, 4), mean)
  names(dim(clim)) <- NULL
  expect_error(
    Ano(var, clim),
    "Provide dimension names on parameter \'var\' and \'clim\' to avoid ambiguity."
  )

  t <- array(rnorm(76), c(1, 3, 4, 3, 2, 2)) 
  names(dim(t)) <- c("mod", "memb", "sdates", "ltime", "lon", "lat")
  c3 <- Clim(t, t, memb = TRUE)$clim_exp # Clim for each member
  c1 <- InsertDim(c3[, 1, ,, ], 1, 1) # clim as if memb=FALSE but identical to member 1              
  names(dim(c1)) <- c("mod",  "ltime", "lon", "lat")
  identical(c1[, , , ], c3[, 1, , , ]) # results in TRUE
  a3 <- Ano(t, c3)  # ano for each member individually                                           
  a1 <- Ano(t, c1)  # ano for first member                                                       

  expect_equal(
    a1[, 1, , , , ],
    a3[, 1, , , , ]
  )

})
