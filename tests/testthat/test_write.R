test_that("basic write sas works", {

  # this x64 file works fine in SAS studio
  tmp <- tempfile(fileext = ".sas7bdat")

  write.sas(mtcars, tmp)
  df <- read.sas(tmp)
  expect_equal(df, mtcars, ignore_attr = TRUE)

  # 32 bit
  write.sas(mtcars, tmp, bit32 = TRUE)
  df <- read.sas(tmp)
  expect_equal(df, mtcars, ignore_attr = TRUE)

})

test_that("basic write sas works", {

  # this x64 file works fine in SAS studio
  data <- iris
  data$Species <- as.character(data$Species)
  tmp <- tempfile(fileext = ".sas7bdat")

  expect_warning(write.sas(data, tmp))
  df <- read.sas(tmp)
  expect_equal(df, data, ignore_attr = TRUE)

  # 32 bit
  expect_warning(write.sas(data, tmp))
  df <- read.sas(tmp)
  expect_equal(df, data, ignore_attr = TRUE)

})
