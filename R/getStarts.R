getStarts <- function(.N,
                      .J,
                      .D,
                      .type = "zeros"
                      ) {

    if (.type == "zeros"){
        starts <- list(alpha = {matrix(0,
                                       nrow = .J,
                                       ncol = 1
                                       )
                            },
                       beta = {matrix(0,
                                      nrow = .J,
                                      ncol = .D
                                      )
                           },
                       x = {matrix(rnorm(.N * .D),
                                   nrow = .N,
                                   ncol = .D
                                   )
                        }
                       )
    } else if (.type == "random") {
        starts <- list(alpha = {matrix(rnorm(.J * 1) * 10,
                                       nrow = .J,
                                       ncol = 1
                                       )
                            },
                       beta = {matrix(rnorm(.J * .D) * 10,
                                      nrow = .J,
                                      ncol = .D
                                      )
                               },
                       x = {matrix(rnorm(.N * .D) * 1,
                                   nrow = .N,
                                   ncol = .D
                                   )
                        }
                       )
    } else {
        stop("Unknown type.")
    }


    return(starts)
}
