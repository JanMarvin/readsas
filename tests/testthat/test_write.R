test_that("basic write sas works", {

  # this x64 file works fine in SAS studio
  tmp <- tempfile(fileext = ".sas7bdat")

  write.sas(mtcars, tmp)
  df <- read.sas(tmp)
  expect_equal(df, mtcars, ignore_attr = TRUE)

  # # 32 bit
  # write.sas(mtcars, tmp, bit32 = TRUE)
  # df <- read.sas(tmp)
  # expect_equal(df, mtcars, ignore_attr = TRUE)

})

test_that("basic write sas works", {

  # this x64 file works fine in SAS studio
  data <- iris
  data$Species <- as.character(data$Species)
  tmp <- tempfile(fileext = ".sas7bdat")

  write.sas(data, tmp)
  df <- read.sas(tmp)
  expect_equal(df, data, ignore_attr = TRUE)

  # # 32 bit
  # write.sas(data, tmp)
  # df <- read.sas(tmp)
  # expect_equal(df, data, ignore_attr = TRUE)

})

test_that("formats work", {

  dd <- data.frame(
    datetime = as.POSIXct("1950-01-02 23:45:43 12:23:45"),
    date = as.Date(Sys.time()),
    num = 1,
    char = "a"
  )

  tmp <- tempfile(fileext = ".sas7bdat")
  write.sas(dd, tmp, bit32 = FALSE)

  attributes(df)$formats


  df <- read.sas(tmp)
  expect_equal(dd, df, ignore_attr = TRUE)

  exp <- c("DATETIME", "DATE", "BEST", "")
  got <- attributes(df)$formats
  expect_equal(exp, got)

})
