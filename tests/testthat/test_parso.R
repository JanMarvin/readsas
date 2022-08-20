test_that("parso external tests", {

  parso_files <- c("all_rand_normal", "all_rand_normal_with_deleted", "all_rand_normal_with_deleted2",
                   "all_rand_uniform", "almost_rand", "comp_deleted","data_page_with_deleted",
                   "doubles", "doubles2", "extend_no", "extend_yes", "file_with_label",
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
    # hv <- haven::read_sas(sas7bdat)
    exp <- read.csv(csv)

    expect_true(all.equal(got, exp, check.attributes = FALSE))
  }

})
