
### basic read
fl <- system.file("extdata", "cars.sas7bdat", package = "readsas")

ds <- read.sas(fl)

test_that("read sas", {
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
  got <- read.sas(fl, select.rows = c(1, 2))
  expect_equal(exp, got, check.attributes = FALSE)

})


### select rows
fl <- system.file("extdata", "mtcars.sas7bdat", package = "readsas")

dd <- read.sas(fl, select.rows = c(2:5), rownames = TRUE)

test_that("select.rows", {
  expect_true(all.equal(dd, mtcars[2:5,], check.attributes = FALSE))
})


### select cols
fl <- system.file("extdata", "mtcars.sas7bdat", package = "readsas")

dd <- read.sas(fl, select.cols = c("VAR1", "mpg", "hp"), rownames = TRUE)

test_that("select.rows", {
  expect_true(all.equal(dd, mtcars[,c("mpg", "hp")], check.attributes = FALSE))
})


### select cols & rows - read rows 2 to 5
fl <- system.file("extdata", "mtcars.sas7bdat", package = "readsas")

dd <- read.sas(fl, select.cols = c("VAR1", "mpg", "hp"), select.rows = c(2:5), rownames = TRUE)

test_that("select.rows", {
  expect_true(all.equal(dd, mtcars[2:5,c("mpg", "hp")], check.attributes = FALSE))
})


### select cols & rows pt 2 - read rows 2 and 5
fl <- system.file("extdata", "mtcars.sas7bdat", package = "readsas")

dd <- read.sas(fl, select.cols = c("VAR1", "mpg", "hp"), select.rows = c(2,5), rownames = TRUE)

test_that("select.rows", {
  expect_true(all.equal(dd, mtcars[c(2,5), c("mpg", "hp")], check.attributes = FALSE))
})


### select cols & rows pt 3 - read zero rows
fl <- system.file("extdata", "mtcars.sas7bdat", package = "readsas")

dd <- read.sas(fl, select.cols = c("VAR1", "mpg", "hp"), select.rows = c(0), rownames = TRUE)

test_that("select.rows", {
  expect_true(all.equal(dd, mtcars[NULL, c("mpg", "hp")], check.attributes = FALSE))
})


test_that("convert time and date", {

  exp <- structure(11033, class = "Date")
  got <- convert_to_date(14686)
  expect_equal(exp, got)

  exp <- structure(1352519359, tzone = "UTC", class = c("POSIXct", "POSIXt"))
  got <- convert_to_datetime(1668138559)
  expect_equal(exp, got)

})

