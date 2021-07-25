
context("Reading datasets")


### basic read
fl <- system.file("extdata", "cars.sas7bdat", package="readsas")

ds <- read.sas(fl)

test_that("read sas", {
  expect_true(all.equal(cars, ds, check.attributes = FALSE))
})


### select rows
fl <- system.file("extdata", "mtcars.sas7bdat", package="readsas")

dd <- read.sas(fl, select.rows = c(2,5))[-1]


test_that("select.rows", {
  expect_true(all.equal(dd, mtcars[2:5,], check.attributes = FALSE))
})


# ### select rows # FIXME: crashing
# fl <- system.file("extdata", "mtcars.sas7bdat", package="readsas")
#
# dd <- read.sas(fl, select.cols = "MPG", debug = F)
#
# test_that("select.rows", {
#   expect_true(all.equal(dd, mtcars[2:5,], check.attributes = FALSE))
# })
