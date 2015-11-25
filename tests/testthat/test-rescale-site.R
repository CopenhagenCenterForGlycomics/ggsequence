context("Rescaling sites")

rescale_site = ggsequence:::rescale_site

test_that("An ungapped sequence doesn't change site", {
    sequence = "MNMNMN"
    expect_equal( rescale_site( sequence , 1), 1 )
    expect_equal( rescale_site( sequence , 3), 3 )
    expect_equal( rescale_site( sequence , 6), 6 )
})

test_that("An single gap in the alignment rescale sites", {
    sequence = "MNM-MN"
    expect_equal( rescale_site( sequence , 1), 1 )
    expect_equal( rescale_site( sequence , 2), 2 )
    expect_equal( rescale_site( sequence , 3), 3 )
    expect_equal( rescale_site( sequence , 4), 5 )
    expect_equal( rescale_site( sequence , 5), 6 )
})

test_that("Multiple gaps in the alignment rescale sites", {
    sequence = "MN-M-MN"
    expect_equal( rescale_site( sequence , 1), 1 )
    expect_equal( rescale_site( sequence , 2), 2 )
    expect_equal( rescale_site( sequence , 3), 4 )
    expect_equal( rescale_site( sequence , 4), 6 )
    expect_equal( rescale_site( sequence , 5), 7 )
})

test_that("Starting gaps rescale sites", {
    sequence = "-MNMMN"
    expect_equal( rescale_site( sequence , 1), 2 )
    expect_equal( rescale_site( sequence , 2), 3 )
    expect_equal( rescale_site( sequence , 3), 4 )
    expect_equal( rescale_site( sequence , 4), 5 )
    expect_equal( rescale_site( sequence , 5), 6 )
})

test_that("Longer starting gaps rescale sites", {
    sequence = "--MNMMN"
    expect_equal( rescale_site( sequence , 1), 3 )
    expect_equal( rescale_site( sequence , 2), 4 )
    expect_equal( rescale_site( sequence , 3), 5 )
    expect_equal( rescale_site( sequence , 4), 6 )
    expect_equal( rescale_site( sequence , 5), 7 )
})

test_that("Ending gaps don't rescale sites", {
    sequence = "MNMMN-"
    expect_equal( rescale_site( sequence , 1), 1 )
    expect_equal( rescale_site( sequence , 2), 2 )
    expect_equal( rescale_site( sequence , 3), 3 )
    expect_equal( rescale_site( sequence , 4), 4 )
    expect_equal( rescale_site( sequence , 5), 5 )
})

test_that("Longer ending gaps don't rescale sites", {
    sequence = "MNMMN--"
    expect_equal( rescale_site( sequence , 1), 1 )
    expect_equal( rescale_site( sequence , 2), 2 )
    expect_equal( rescale_site( sequence , 3), 3 )
    expect_equal( rescale_site( sequence , 4), 4 )
    expect_equal( rescale_site( sequence , 5), 5 )
})

test_that("A combination of gaps rescale sites", {
    sequence = "-M-N-M-M-N-"
    expect_equal( rescale_site( sequence , 1), 2 )
    expect_equal( rescale_site( sequence , 2), 4 )
    expect_equal( rescale_site( sequence , 3), 6 )
    expect_equal( rescale_site( sequence , 4), 8 )
    expect_equal( rescale_site( sequence , 5), 10 )
})

