test_that("parso external tests", {

  # # disable for this branch
  # parso_files <- c(
  #   "all_rand_normal", "all_rand_normal_with_deleted", "all_rand_normal_with_deleted2",
  #   "all_rand_uniform", "almost_rand", "comp_deleted", "data_page_with_deleted",
  #   "doubles", "doubles2", "extend_no", "extend_yes", "file_with_label",
  #   "fileserrors", "int_only", "int_only_partmissing",
  #   "mix_and_missing",                                # character "" is not imported as NA
  #   "mix_data_misc", "mix_data_with_longchar",
  #   "mix_data_with_missing_char",                     # character "" is not imported as NA
  #   "mix_data_with_partmissing_char",                 # character "" is not imported as NA
  #   "mixed_data_one", "mixed_data_two",
  #   "only_datetime",
  #   # "test-columnar",                                  # formatted pct date, datetime and time
  #   # "time_formats",                                   # formatted time
  #   "tmp868_14"
  # )
  #
  # for (fl in parso_files) {
  #
  #   # message(fl)
  #   sas7bdat <- sprintf("https://github.com/epam/parso/raw/master/src/test/resources/sas7bdat/%s.sas7bdat", fl)
  #   csv      <- sprintf("https://github.com/epam/parso/raw/master/src/test/resources/csv/%s.csv", fl)
  #   # # testing with local parso
  #   # sas7bdat <- sprintf("~/Source/parso/src/test/resources/sas7bdat/%s.sas7bdat", fl)
  #   # csv      <- sprintf("~/Source/parso/src/test/resources/csv/%s.csv", fl)
  #
  #   empty_to_na <- FALSE
  #
  #   if (fl %in% c("mix_and_missing", "mix_data_with_missing_char", "mix_data_with_partmissing_char"))
  #     empty_to_na <- TRUE
  #
  #   got <- read.sas(sas7bdat, debug = FALSE, empty_to_na = empty_to_na)
  #   exp <- read.csv(csv, na.strings = "NA")
  #
  #   # file is written with character missings. These are converted to integers. read.sas otherwise
  #   # imports characters and read.csv imports integers and logical for entirely missing columns. This
  #   # simply fixes the comparison, otherwise warnings regarding NA conversion would appear.
  #   if (fl %in% c("mix_and_missing", "mix_data_with_missing_char", "mix_data_with_partmissing_char")) {
  #     got <- as.data.frame(apply(got, 2, as.integer))
  #     exp <- as.data.frame(apply(exp, 2, as.integer))
  #   }
  #
  #   expect_true(all.equal(got, exp, check.attributes = FALSE))
  # }

})

test_that("parso external char tests", {

  parso_files <- c(#"charset_aara",
    #"charset_acro",
    "charset_acyr", "charset_agrk",
    # "charset_aheb",
    # "charset_aice",
    "charset_ansi", "charset_arab",
    # "charset_arma",
    "charset_arom",
    # "charset_atha",
    "charset_atur",
    "charset_aukr", "charset_big5", "charset_cyrl", "charset_gbke",
    "charset_grek", "charset_hebr", "charset_j932", "charset_j942",
    "charset_jeuc", "charset_jiso", "charset_keuc", "charset_kiso",
    "charset_kpce", "charset_kwin", "charset_lat1", "charset_lat2",
    "charset_lat3", "charset_lat4", "charset_lat5", "charset_lat7",
    "charset_lat9", "charset_p852", "charset_p857", "charset_p858",
    "charset_p860", "charset_p862", "charset_p864", "charset_p874",
    "charset_sjis",
    # "charset_sjs4",
    "charset_thai", "charset_utf8",
    "charset_wara", "charset_wbal", "charset_wcyr", "charset_wgrk",
    "charset_wheb", "charset_wlt1", "charset_wlt2", "charset_wtur",
    "charset_wvie", "charset_yeuc", "charset_ywin", "charset_zeuc"#,
    # "charset_zpce",
    # "charset_zwin"
  )

  for (fl in parso_files) {

    # message(fl)
    sas7bdat <- sprintf("https://github.com/epam/parso/raw/master/src/test/resources/sas7bdat/%s.sas7bdat", fl)
    csv      <- sprintf("https://github.com/epam/parso/raw/master/src/test/resources/csv/%s.csv", fl)

    got <- read.sas(sas7bdat, debug = FALSE)
    exp <- read.csv(csv, na.strings = "NA")

    expect_true(all.equal(got, exp, check.attributes = FALSE))
  }

})

test_that("test external sas7bdat", {

  sas7bdat_files <- c(
    "charactervalues",
    "datetimevalues",
    "datevalues",
    "floatvalues",
    "intvalues",
    "mixedvalues",
    "mixedvalues_compressed_binary",
    "mixedvalues_compressed_char",
    "mixedvalues_compressed_yes",
    "mixedvalues_empty",
    "specialcharactervalues",
    "timevalues"
  )

  for (fl in sas7bdat_files) {

    # message(fl)
    sas7bdat <- sprintf("https://bitbucket.org/jaredhobbs/sas7bdat/raw/18cbd14407918c1aa90f136c1d6c5d83f307dba0/tests/data/%s.sas7bdat", fl)

    exp <- "SAS FILE"
    got <- read.sas(sas7bdat, debug = FALSE)
    got <- attributes(got)$sasfile
    expect_equal(exp, got)
  }

})
