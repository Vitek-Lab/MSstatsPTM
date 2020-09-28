test_that("annotSite() annotates PTM sites", {
    expect_equal(annotSite(10, "K"), "K10")
    expect_equal(annotSite(10, "K", 3L), "K010")
})
