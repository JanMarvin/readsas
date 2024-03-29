test_that("parso external tests", {

  parso_files <- c(
    "all_rand_normal", "all_rand_normal_with_deleted", "all_rand_normal_with_deleted2",
    "all_rand_uniform", "almost_rand", "comp_deleted", "data_page_with_deleted",
    "doubles", "doubles2", "extend_no", "extend_yes", "file_with_label",
    "fileserrors", "int_only", "int_only_partmissing",
    "mix_and_missing",                                # character "" is not imported as NA
    "mix_data_misc", "mix_data_with_longchar",
    "mix_data_with_missing_char",                     # character "" is not imported as NA
    "mix_data_with_partmissing_char",                 # character "" is not imported as NA
    "mixed_data_one", "mixed_data_two",
    "only_datetime",
    # "test-columnar",                                  # formatted pct date, datetime and time
    # "time_formats",                                   # formatted time
    "tmp868_14"
  )

  for (fl in parso_files) {

    message(fl)
    sas7bdat <- sprintf("https://github.com/epam/parso/raw/master/src/test/resources/sas7bdat/%s.sas7bdat", fl)
    csv      <- sprintf("https://github.com/epam/parso/raw/master/src/test/resources/csv/%s.csv", fl)
    # # testing with local parso
    # sas7bdat <- sprintf("~/Source/parso/src/test/resources/sas7bdat/%s.sas7bdat", fl)
    # csv      <- sprintf("~/Source/parso/src/test/resources/csv/%s.csv", fl)

    empty_to_na <- FALSE

    if (fl %in% c("mix_and_missing", "mix_data_with_missing_char", "mix_data_with_partmissing_char"))
      empty_to_na <- TRUE

    got <- read.sas(sas7bdat, debug = FALSE, empty_to_na = empty_to_na)
    exp <- read.csv(csv, na.strings = "NA")

    # file is written with character missings. These are converted to integers. read.sas otherwise
    # imports characters and read.csv imports integers and logical for entirely missing columns. This
    # simply fixes the comparison, otherwise warnings regarding NA conversion would appear.
    if (fl %in% c("mix_and_missing", "mix_data_with_missing_char", "mix_data_with_partmissing_char")) {
      got <- as.data.frame(apply(got, 2, as.integer))
      exp <- as.data.frame(apply(exp, 2, as.integer))
    }

    expect_true(all.equal(got, exp, check.attributes = FALSE))
  }

})
