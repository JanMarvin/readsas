
context("Reading datasets")


fl <- system.file("extdata", "cars.sas7bdat", package="readsas")

ds <- read.sas(fl)

test_that("compare", {
  expect_true(all.equal(cars, ds, check.attributes = FALSE))
})
