library(testthat)
library(CVE)
library(pingr)

test_that("oncotator example file is in the right format",
          {
            expect_that(oncotator_example, is.data.frame)
          })

test_that("ping local server for Shiny app (less than 1 sec)",
          {
          expect_less_than(ping("127.0.0.1",count = 1), 1)
          })


