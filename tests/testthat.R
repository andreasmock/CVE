library(testthat)
library(CVE)
library(pingr)

test_that("oncotator example file is in the right format",
          {
            expect_that(oncotator_example, is.data.frame)
          })

test_that("ping local server for Shiny app",
          {
          expect_equal(sum(ping("127.0.0.1")>1), 0)
          })


