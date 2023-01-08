test_that("basic write sas works", {

  # this file works fine in SAS studio
  tmp <- tempfile(fileext = ".sas7bdat")
  write.sas(mtcars, tmp)

  df <- read.sas(tmp)

  expect_equal(df, mtcars, ignore_attr = TRUE)

})

test_that("basic write sas works", {

  # this file works fine in SAS studio
  data <- iris
  data$Species <- as.character(data$Species)

  tmp <- tempfile(fileext = ".sas7bdat")
  expect_warning(write.sas(data, tmp))

  df <- read.sas(tmp)

  expect_equal(df, data, ignore_attr = TRUE)

})
