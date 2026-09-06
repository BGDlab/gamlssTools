# `compare_centile_fans()` builds its plot by transplanting the second model's line
# (and peak) layers onto the first model's plot, then stripping the color mappings
# make_centile_fan() put there so that color can carry the MODEL instead and facet_var
# can become facets. The tests below check the structure that surgery produces -- one
# line layer per model, distinct fixed colors, no leftover mapped color, one panel per
# level -- rather than just that the call returned without error. The failure modes that
# motivate them (a plot showing one model twice, both models in the same color, or a
# leftover facet_var mapping fighting the model scale) all render without any error.

skip_if_not_installed("gamlss")

cmp_fixture <- function(seed = 11, n = 150) {
  set.seed(seed)
  d <- data.frame(Age = runif(n, 1, 20),
                  Sex = factor(sample(c("M", "F"), n, TRUE)))
  d$Pheno <- 5 + log(d$Age) + (d$Sex == "M") * 0.4 + rnorm(n, 0, 0.4)
  ctl <- gamlss::gamlss.control(trace = FALSE)
  list(d  = d,
       m1 = gamlss::gamlss(Pheno ~ Age + Sex, data = d, control = ctl),
       m2 = gamlss::gamlss(Pheno ~ pb(Age) + Sex, data = d, control = ctl))
}

line_layers  <- function(p) Filter(function(l) inherits(l$geom, "GeomLine"), p$layers)
point_layers <- function(p) Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)
fan_layers   <- function(p) Filter(function(l) !is.null(l$data) && is.data.frame(l$data) &&
                                     nrow(l$data) > 2, line_layers(p))
quiet <- function(expr) suppressMessages(utils::capture.output(res <- expr))

test_that("both models' fans are drawn solid, in their own color, at the given alpha", {
  f <- cmp_fixture()
  quiet(p <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                  desiredCentiles = c(0.05, 0.5, 0.95)))
  expect_s3_class(p, "ggplot")

  #two fans plus the invisible legend-carrying layer
  lines <- line_layers(p)
  expect_length(lines, 3)

  fans <- fan_layers(p)
  expect_length(fans, 2)
  expect_setequal(vapply(fans, function(l) l$aes_params$colour, character(1)),
                  c("blue", "red"))
  expect_true(all(vapply(fans, function(l) l$aes_params$linetype, character(1)) == "solid"))
  expect_true(all(vapply(fans, function(l) l$aes_params$alpha, numeric(1)) == 0.6))

  #the two fans hold DIFFERENT predictions: a silently duplicated model would tie
  expect_false(isTRUE(all.equal(fans[[1]]$data, fans[[2]]$data)))
})

test_that("alpha and model_colors are honoured", {
  f <- cmp_fixture()
  quiet(p <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                  model_colors = c("#1B9E77", "#D95F02"), alpha = 0.35))
  fans <- fan_layers(p)
  expect_setequal(vapply(fans, function(l) l$aes_params$colour, character(1)),
                  c("#1B9E77", "#D95F02"))
  expect_true(all(vapply(fans, function(l) l$aes_params$alpha, numeric(1)) == 0.35))
})

test_that("no layer is left mapping color, so the model scale has no rival", {
  f <- cmp_fixture()
  quiet(p <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex", get_peaks = TRUE))

  #every layer but the legend layer must carry a FIXED color, not a mapped one:
  #a leftover facet_var mapping would need a second discrete scale on the same aesthetic
  real <- utils::head(p$layers, -1)
  expect_true(all(vapply(real, function(l) is.null(l$mapping$colour), logical(1))))
  expect_true(all(vapply(real, function(l) is.null(l$mapping$fill), logical(1))))

  legend_layer <- p$layers[[length(p$layers)]]
  expect_false(is.null(legend_layer$mapping$colour))
})

test_that("the model legend carries one key per model, named after the objects", {
  f <- cmp_fixture()
  quiet(p <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex"))

  legend_data <- p$layers[[length(p$layers)]]$data
  expect_equal(levels(legend_data$model), c("f$m1", "f$m2"))
  #the legend layer must draw nothing itself
  expect_true(all(is.na(legend_data$x)))

  quiet(named <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                      model_names = c("linear", "spline")))
  expect_equal(levels(named$layers[[length(named$layers)]]$data$model),
               c("linear", "spline"))

  #the same model passed twice still needs two distinct keys
  quiet(dup <- compare_centile_fans(f$m1, f$m1, f$d, "Age", "Sex"))
  expect_equal(levels(dup$layers[[length(dup$layers)]]$data$model),
               c("f$m1 (1)", "f$m1 (2)"))
})

test_that("facet_var becomes facets, and only when there is something to facet", {
  f <- cmp_fixture()

  quiet(p <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex"))
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(names(p$facet$params$facets), "Sex")
  #the panels must actually be built, one per level
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), nlevels(f$d$Sex))

  #averaging leaves one fan per model and nothing to facet by
  quiet(avg <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex", average_over = TRUE))
  expect_s3_class(avg$facet, "FacetNull")

  #and no facet_var means no facets either
  quiet(none <- compare_centile_fans(f$m1, f$m2, f$d, "Age"))
  expect_s3_class(none$facet, "FacetNull")

  quiet(free <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex", facet_scales = "free_y"))
  expect_true(free$facet$params$free$y)
})

test_that("points are drawn once, from the first model's data, in a neutral color", {
  f <- cmp_fixture()
  quiet(p <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                  show_points = TRUE, get_peaks = FALSE,
                                  point_color = "grey30",
                                  label_centiles = "label"))

  #exactly one geom_point layer, holding every row of the data
  pts <- point_layers(p)
  expect_length(pts, 1)
  expect_equal(nrow(pts[[1]]$data), nrow(f$d))
  expect_equal(pts[[1]]$aes_params$colour, "grey30")

  #labels come from the first model's fan only, so one label per centile per level
  #(three centiles by default, not make_centile_fan's nine)
  labels <- Filter(function(l) inherits(l$geom, "GeomTextRepel"), p$layers)
  expect_length(labels, 1)
  expect_equal(nrow(labels[[1]]$data), 3 * nlevels(f$d$Sex))
})

test_that("peaks are marked with a different shape and each model's color", {
  f <- cmp_fixture()
  quiet(p <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                  show_points = FALSE, get_peaks = TRUE))
  peaks <- point_layers(p)
  expect_setequal(vapply(peaks, function(l) as.numeric(l$aes_params$shape), numeric(1)),
                  c(18, 17))
  expect_setequal(vapply(peaks, function(l) l$aes_params$colour, character(1)),
                  c("blue", "red"))
})

test_that("a shared sim grid gives the same fans as letting each call build its own", {
  f <- cmp_fixture()
  quiet(grid <- sim_grid(f$d, "Age", "Sex", f$m1))
  quiet(own <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex"))
  fan_data <- function(p) lapply(fan_layers(p), function(l) l$data)

  quiet(both <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                     sim_grid_list = grid, sim_grid_list2 = grid))
  expect_equal(fan_data(both), fan_data(own))

  #one grid supplied is reused for the second model, matching how df2 falls back to df
  quiet(one <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex", sim_grid_list = grid))
  expect_equal(fan_data(one), fan_data(own))
})

test_that("arguments set per model are refused when passed through `...`", {
  f <- cmp_fixture()
  expect_error(compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                    sim_data_list = list()),
               "sim_grid_list")
  expect_error(compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                    color_manual = "red"),
               "model_colors")
  expect_error(compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                    point_color_manual = "red"),
               "point_color")
})

test_that("the percentile legend's keys are neutral, and only when it is drawn", {
  f <- cmp_fixture()

  #without this the keys inherit the first model's color, reading as though line
  #thickness belonged to that model
  quiet(leg <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex",
                                    label_centiles = "legend"))
  expect_true("linewidth" %in% names(leg$guides$guides))
  expect_equal(leg$guides$guides$linewidth$params$override.aes$colour, "black")

  #a plot-level guides() overrides the scale's own guide = "none", so setting it for the
  #other options would bring the percentile legend back for plots that asked not to have it
  for (lc in c("label", "none")) {
    quiet(p <- compare_centile_fans(f$m1, f$m2, f$d, "Age", "Sex", label_centiles = lc))
    expect_false("linewidth" %in% names(p$guides$guides), info = lc)
  }
})
