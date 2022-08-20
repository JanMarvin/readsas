
fl <- system.file("extdata", "cars.sas7bdat", package = "readsas")

ds <- read.sas(fl)

test_that("compare", {
  expect_true(all.equal(cars, ds, check.attributes = FALSE))
})


test_that("read file with deleted rows", {

  # read full file
  fl <- system.file("extdata", "test.sas7bdat", package = "readsas")
  exp <- data.frame(x = 1:3)
  got <- read.sas(fl)
  expect_equal(exp, got, check.attributes = FALSE)

  # read file with deleted: delete x = 2
  fl <- system.file("extdata", "test2.sas7bdat", package = "readsas")
  exp <- data.frame(x = 1:3)[c(1, 3), , drop = FALSE]
  got <- read.sas(fl)
  expect_equal(exp, got, check.attributes = FALSE)

  # read file with deleted: delete x > 1
  fl <- system.file("extdata", "test3.sas7bdat", package = "readsas")
  exp <- data.frame(x = 1:3)[c(1), , drop = FALSE]
  got <- read.sas(fl)
  expect_equal(exp, got, check.attributes = FALSE)

  # read file with deleted: delete x = 1
  fl <- system.file("extdata", "test4.sas7bdat", package = "readsas")
  exp <- data.frame(x = 1:3)[-c(1), , drop = FALSE]
  got <- read.sas(fl)
  expect_equal(exp, got, check.attributes = FALSE)

  # read file with deleted: delete x < 3
  fl <- system.file("extdata", "test5.sas7bdat", package = "readsas")
  exp <- data.frame(x = 1:3)[-c(1, 2), , drop = FALSE]
  got <- read.sas(fl)
  expect_equal(exp, got, check.attributes = FALSE)

  # skip row
  fl <- system.file("extdata", "test2.sas7bdat", package = "readsas")
  exp <- data.frame(x = 1:3)[c(1), , drop = FALSE]
  expect_warning(
    got <- read.sas(fl, rowcount = 2),
    "User requested to read 2 of 3 rows"
  )
  expect_equal(exp, got, check.attributes = FALSE)

})
