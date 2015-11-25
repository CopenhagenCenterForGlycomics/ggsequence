context("Rescaling sites")

rescale_site = ggsequence:::rescale_site

test_that("An ungapped sequence doesn't change site", {
    sequence = "MNMNMN"
    expect_equal( rescale_site( sequence , 1), 1 )
    expect_equal( rescale_site( sequence , 3), 3 )
    expect_equal( rescale_site( sequence , 6), 6 )
})

test_that("An single gap in the alignment", {
    sequence = "MNM-MN"
    expect_equal( rescale_site( sequence , 1), 1 )
    expect_equal( rescale_site( sequence , 2), 2 )
    expect_equal( rescale_site( sequence , 3), 3 )
    expect_equal( rescale_site( sequence , 4), 5 )
    expect_equal( rescale_site( sequence , 5), 6 )
})