
context("Reading datasets")


fl <- system.file("extdata", "cars.sas7bdat", package = "readsas")

ds <- read.sas(fl)

test_that("compare", {
  expect_true(all.equal(cars, ds, check.attributes = FALSE))
})

test_that("parso external tests", {

  parso_files <- c("all_rand_normal", "all_rand_normal_with_deleted", "all_rand_normal_with_deleted2",
                   "all_rand_uniform", "almost_rand",
                   "comp_deleted",                                     # deleted row at end of file
                   "data_page_with_deleted",
                   "doubles", "doubles2", "extend_no",
                   # "extend_yes",                                     # varnames applied incorrectly
                   "file_with_label",
                   "fileserrors", "int_only", "int_only_partmissing",
                   # "mix_and_missing",                                # character "" is not imported as NA
                   "mix_data_misc", "mix_data_with_longchar",
                   # "mix_data_with_missing_char",                     # character "" is not imported as NA
                   # "mix_data_with_partmissing_char",                 # character "" is not imported as NA
                   "mixed_data_one", "mixed_data_two",
                   "only_datetime",
                   # "test-columnar",                                  # formatted pct date, datetime and time
                   # "time_formats",                                   # formatted time
                   "tmp868_14")



  for (fl in parso_files) {

    message(fl)
    sas7bdat <- sprintf("https://github.com/epam/parso/raw/master/src/test/resources/sas7bdat/%s.sas7bdat", fl)
    csv      <- sprintf("https://github.com/epam/parso/raw/master/src/test/resources/csv/%s.csv", fl)

    got <- read.sas(sas7bdat)
    # got <- haven::read_sas(sas7bdat)
    exp <- read.csv(csv)

    expect_true(all.equal(got, exp, check.attributes = FALSE))
  }

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
