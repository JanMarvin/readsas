
context("Reading datasets")


fl <- system.file("extdata", "cars.sas7bdat", package="readsas")

ds <- read.sas(fl)
ds1 <- read.sas(fl, select.rows = 4)
ds2 <- read.sas(fl, select.rows = c(2,4))

test_that("compare", {
  expect_true(all.equal(cars, ds, check.attributes = FALSE))
  expect_true(all.equal(cars[1:4,], ds1, check.attributes = FALSE))
  expect_true(all.equal(cars[2:4,], ds2, check.attributes = FALSE))
})
