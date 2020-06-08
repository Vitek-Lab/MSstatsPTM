test_that("annot_site() annotates PTM sites", {
    expect_equal(annot_site(10, "K"), "K10")
    expect_equal(annot_site(10, "K", 3L), "K010")
})
