# =========================================================
# Shiny App: t, F, Chi-square, Spearman and Mann-Whitney U
# distributions derived by simulation
# (with controlled resampling + previous sample overlay)
# =========================================================

# https://hbctraining.github.io/Training-modules/RShiny/lessons/shinylive.html
# Run the shinylive::export line to populate the docs folder
# so that shinylive works from github
#shinylive::export(appdir = "../DerivingClassicDistributions/", destdir = "docs")
#httpuv::runStaticServer("docs/", port = 8008)

# ---------------------------------------------------------
# Safari-safe plotting
# Draws a base-R plot as a PNG image. Stand-in for renderPlot(), which
# currently fails in Safari under Shinylive. Use with imageOutput().
# ---------------------------------------------------------
render_base <- function(id, session, draw) {
  renderImage({
    cd <- session$clientData
    w  <- cd[[paste0("output_", id, "_width")]]
    h  <- cd[[paste0("output_", id, "_height")]]
    req(w > 0, h > 0)
    pr <- cd$pixelratio
    if (is.null(pr)) pr <- 1
    f <- tempfile(fileext = ".png")
    png(f, width = w * pr, height = h * pr, res = 72 * pr)
    tryCatch(draw(), finally = dev.off())
    list(src = f, contentType = "image/png", width = w, height = h)
  }, deleteFile = TRUE)
}

# ---------------------------------------------------------
# Histogram of the current sample, with the previous sample
# drawn underneath in grey (used by the t, chi-square and Spearman plots)
# ---------------------------------------------------------
hist_with_prev <- function(vals, prev_vals, breaks, main, xlab, ylim_extra = 0) {
  h_curr <- hist(vals, breaks = breaks, plot = FALSE)
  ylim   <- c(0, max(c(h_curr$density, ylim_extra)) * 1.15)
  if (!is.null(prev_vals)) {
    hist(prev_vals, probability = TRUE, breaks = breaks,
         col = rgb(0.7, 0.7, 0.7, 0.4), border = NA,
         main = main, xlab = xlab, ylim = ylim)
    hist(vals, probability = TRUE, breaks = breaks,
         col = "#569BBD", border = "black", add = TRUE)
  } else {
    hist(vals, probability = TRUE, breaks = breaks,
         col = "#569BBD", border = "black",
         main = main, xlab = xlab, ylim = ylim)
  }
}

# Theoretical (null) density of a correlation coefficient with sample size n.
# Exact for Pearson's r; a close approximation for Spearman's r_S.
dcor_null <- function(r, n) {
  (1 - r^2)^((n - 4) / 2) / beta(1/2, (n - 2) / 2)
}

# Critical value of a correlation coefficient (alpha = 0.05, two-tailed),
# found from the t distribution with n - 2 df
r_crit <- function(n) {
  tc <- qt(0.975, df = n - 2)
  tc / sqrt(tc^2 + n - 2)
}

# Dashed line at a critical value (alpha = 0.05), labelled with its value.
# side = "right" puts the label to the right of the line, "left" to the left.
add_crit <- function(x, label, side = "right") {
  abline(v = x, lty = 2, lwd = 1.5, col = "grey25")
  text(x, par("usr")[4] * 0.93, label, pos = if (side == "right") 4 else 2,
       cex = 0.9, col = "grey20", font = 2)
}

# Dotted orange line at the critical value found from the simulated values
# themselves (empirical), labelled just below the theoretical label.
emp_col <- "#E07B39"
add_emp <- function(x, label, side = "right") {
  abline(v = x, lty = 3, lwd = 2.5, col = emp_col)
  text(x, par("usr")[4] * 0.83, label, pos = if (side == "right") 4 else 2,
       cex = 0.9, col = emp_col, font = 2)
}

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
    #resample {
      background-color: #569BBD;
      color: white;
      border: none;
    }
    #resample:hover {
      background-color: #3E7C99;
      color: white;
    }
  "))
  ),
  
  titlePanel("Deriving Classical Distributions from Normal Sampling"),
  
  sidebarLayout(
    sidebarPanel(
      
      # -------------------------
      # Distribution selector
      # -------------------------
      selectInput(
        "dist",
        "Choose distribution:",
        choices = c("t distribution" = "t",
                    "F distribution" = "f",
                    "Chi-square distribution" = "chi",
                    "Binomial distribution" = "binom",
                    "Pearson's correlation (r)" = "pear",
                    "Spearman's rank correlation" = "spear",
                    "Mann-Whitney U" = "mwu")
      ),
      
      # -------------------------
      # Resample button
      # -------------------------
      actionButton("resample", "Resample"),
      
      hr(),
      
      # -------------------------
      # Dynamic controls
      # -------------------------
      uiOutput("controls"),
      
      helpText("_______________________________"),
      helpText("Glenn Tattersall, PhD"),
      helpText("For use in BIOL 3P96 - Biostatistics"),
      br()
    ),
    
    
    
    mainPanel(
      imageOutput("distPlot", height = "400px"),
      htmlOutput("explain_text"),
      verbatimTextOutput("explain_code")
    )
  )
)

server <- function(input, output, session) {
  
  # -----------------------------
  # Dynamic UI controls
  # -----------------------------
  output$controls <- renderUI({
    
    if (input$dist == "t") {
      tagList(
        h4("t Distribution Settings"),
        sliderInput("n_t", "Sample size per t experiment (n):",
                    min = 5, max = 100, value = 30),
        sliderInput("n_exp_t", "Number of simulations:",
                    min = 100, max = 5000, value = 1000, step = 100)
      )
      
    } else if (input$dist == "f") {
      tagList(
        h4("F Distribution Settings"),
        sliderInput("n1_f", "Sample size 1 (numerator df = n1-1):",
                    min = 2, max = 100, value = 10),
        sliderInput("n2_f", "Sample size 2 (denominator df = n2-1):",
                    min = 2, max = 100, value = 10),
        sliderInput("n_exp_f", "Number of simulations:",
                    min = 100, max = 5000, value = 1000, step = 100)
      )
      
    } else if (input$dist == "chi") {
      tagList(
        h4("Chi-square Distribution Settings"),
        sliderInput("k_chi", "Degrees of freedom (k):",
                    min = 1, max = 50, value = 5),
        sliderInput("n_exp_chi", "Number of simulations:",
                    min = 100, max = 5000, value = 1000, step = 100)
      )
      
    } else if (input$dist == "binom") {
      tagList(
        h4("Binomial Settings"),
        sliderInput("n_b", "Number of trials (n):",
                    min = 5, max = 100, value = 20),
        sliderInput("p_b", "Probability of success in each trial (p):",
                    min = 0.05, max = 0.95, value = 0.5, step = 0.05),
        sliderInput("n_exp_b", "Number of simulations:",
                    min = 100, max = 5000, value = 1000, step = 100)
      )
      
    } else if (input$dist == "pear") {
      tagList(
        h4("Pearson Correlation Settings"),
        sliderInput("n_r", "Number of pairs of observations (n):",
                    min = 5, max = 100, value = 15),
        sliderInput("n_exp_r", "Number of simulations:",
                    min = 100, max = 5000, value = 1000, step = 100)
      )
      
    } else if (input$dist == "spear") {
      tagList(
        h4("Spearman Correlation Settings"),
        sliderInput("n_sp", "Number of pairs of observations (n):",
                    min = 5, max = 100, value = 15),
        sliderInput("n_exp_sp", "Number of simulations:",
                    min = 100, max = 5000, value = 1000, step = 100)
      )
      
    } else {
      tagList(
        h4("Mann-Whitney U Settings"),
        sliderInput("n1_u", "Sample size of group 1 (n1):",
                    min = 3, max = 20, value = 8),
        sliderInput("n2_u", "Sample size of group 2 (n2):",
                    min = 3, max = 20, value = 8),
        sliderInput("n_exp_u", "Number of simulations:",
                    min = 100, max = 5000, value = 1000, step = 100)
      )
    }
  })
  
  # -----------------------------
  # Store previous samples
  # -----------------------------
  prev_t   <- reactiveVal(NULL)
  prev_f   <- reactiveVal(NULL)
  prev_chi <- reactiveVal(NULL)
  prev_b   <- reactiveVal(NULL)
  prev_r   <- reactiveVal(NULL)
  prev_sp  <- reactiveVal(NULL)
  prev_u   <- reactiveVal(NULL)
  
  # =========================================================
  # Controlled simulation (only updates on Resample click)
  # =========================================================
  
  t_data <- reactive({
    input$resample
    req(input$n_t, input$n_exp_t)
    
    old <- isolate(prev_t())
    
    new <- replicate(input$n_exp_t, {
      x <- rnorm(input$n_t)
      mean(x) / (sd(x) / sqrt(input$n_t))
    })
    
    prev_t(new)
    attr(new, "prev") <- old
    new
  })
  
  
  f_data <- reactive({
    input$resample
    req(input$n1_f, input$n2_f, input$n_exp_f)
    
    old <- isolate(prev_f())
    
    new <- replicate(input$n_exp_f, {
      x1 <- rnorm(input$n1_f)
      x2 <- rnorm(input$n2_f)
      var(x1) / var(x2)
    })
    
    prev_f(new)
    attr(new, "prev") <- old
    new
  })
  
  
  chi_data <- reactive({
    input$resample
    req(input$k_chi, input$n_exp_chi)
    
    old <- isolate(prev_chi())
    
    new <- replicate(input$n_exp_chi, {
      z <- rnorm(input$k_chi)
      sum(z^2)
    })
    
    prev_chi(new)
    attr(new, "prev") <- old
    new
  })
  
  
  b_data <- reactive({
    input$resample
    req(input$n_b, input$p_b, input$n_exp_b)
    
    old <- isolate(prev_b())
    n   <- input$n_b
    p   <- input$p_b
    
    # Each trial is a success with probability p; count the successes
    new <- replicate(input$n_exp_b, sum(runif(n) < p))
    attr(new, "n") <- c(n, p)                   # remember the settings used
    
    prev_b(new)
    attr(new, "prev") <- old
    new
  })
  
  
  r_data <- reactive({
    input$resample
    req(input$n_r, input$n_exp_r)
    
    old <- isolate(prev_r())
    
    # X and Y are unrelated (null hypothesis true)
    new <- replicate(input$n_exp_r, {
      x <- rnorm(input$n_r)
      y <- rnorm(input$n_r)
      cor(x, y)
    })
    
    prev_r(new)
    attr(new, "prev") <- old
    new
  })
  
  
  sp_data <- reactive({
    input$resample
    req(input$n_sp, input$n_exp_sp)
    
    old <- isolate(prev_sp())
    
    # X and Y are unrelated (null hypothesis true)
    new <- replicate(input$n_exp_sp, {
      x <- rank(rnorm(input$n_sp))
      y <- rank(rnorm(input$n_sp))
      cor(x, y)                    # correlation of the ranks = Spearman's r_S
    })
    
    prev_sp(new)
    attr(new, "prev") <- old
    new
  })
  
  
  u_data <- reactive({
    input$resample
    req(input$n1_u, input$n2_u, input$n_exp_u)
    
    old <- isolate(prev_u())
    n1  <- input$n1_u
    n2  <- input$n2_u
    
    # Both groups from the same population (null hypothesis true)
    new <- replicate(input$n_exp_u, {
      y1 <- rnorm(n1)
      y2 <- rnorm(n2)
      R1 <- sum(rank(c(y1, y2))[1:n1])          # rank sum of group 1
      U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
      U2 <- n1 * n2 - U1
      max(U1, U2)                               # U = the larger of U1 and U2
    })
    attr(new, "n") <- c(n1, n2)                 # remember the sample sizes used
    
    prev_u(new)
    attr(new, "prev") <- old
    new
  })
  
  
  # -----------------------------
  # Main plot
  # -----------------------------
  output$distPlot <- render_base("distPlot", session, function() {
    
    if (input$dist == "t") {
      
      t_vals    <- t_data()
      prev_vals <- attr(t_vals, "prev")
      df        <- input$n_t - 1
      
      breaks <- pretty(range(c(t_vals, prev_vals), na.rm = TRUE), n = 30)
      hist_with_prev(t_vals, prev_vals, breaks,
                     main = paste("t distribution (n =", input$n_t, ", df =", df, ")"),
                     xlab = "t statistic")
      curve(dt(x, df = df), add = TRUE, col = "red", lwd = 2)
      tc <- qt(0.975, df = df)
      add_crit(-tc, sprintf("t = %.2f", -tc), "left")
      add_crit( tc, sprintf("t = %.2f",  tc), "right")
      te <- quantile(t_vals, c(0.025, 0.975))
      add_emp(te[1], sprintf("sim. %.2f", te[1]), "left")
      add_emp(te[2], sprintf("sim. %.2f", te[2]), "right")
      
    } else if (input$dist == "f") {
      
      f_vals    <- f_data()
      prev_vals <- attr(f_vals, "prev")
      df1       <- input$n1_f - 1
      df2       <- input$n2_f - 1
      
      # Fixed 0-10 axis with bins of width 0.25. Very large F values (common
      # when sample size 2 is small) are left off the axis instead of
      # stretching the bins, and the bars are scaled by ALL simulations so
      # they still line up with the red theoretical curve.
      xmax   <- 10
      breaks <- seq(0, xmax, by = 0.25)
      bin_density <- function(v) {
        h <- hist(v[v <= xmax], breaks = breaks, plot = FALSE)
        h$density <- h$counts / (length(v) * diff(breaks)[1])
        h
      }
      h_curr <- bin_density(f_vals)
      ylim   <- c(0, max(c(h_curr$density, 0.05)) * 1.15)   # 0.05 floor avoids a zero-height axis
      
      title_txt <- paste("F distribution (df1 =", df1, ", df2 =", df2, ")")
      
      plot(NULL, xlim = c(0, xmax), ylim = ylim,
           main = title_txt, xlab = "F statistic", ylab = "Density")
      if (!is.null(prev_vals)) {
        plot(bin_density(prev_vals), freq = FALSE, add = TRUE,
             col = rgb(0.7, 0.7, 0.7, 0.4), border = NA)
      }
      plot(h_curr, freq = FALSE, add = TRUE,
           col = "#569BBD", border = "black")
      
      # Say how many values fell beyond the axis
      n_off <- sum(f_vals > xmax)
      if (n_off > 0) {
        mtext(paste0(n_off, " of ", length(f_vals), " values (",
                     round(100 * n_off / length(f_vals), 1),
                     "%) are above ", xmax, " and not shown"),
              side = 3, line = 0.2, cex = 0.85, col = "grey40")
      }
      
      curve(df(x, df1 = df1, df2 = df2), add = TRUE, col = "red", lwd = 2)
      fc <- qf(0.95, df1 = df1, df2 = df2)
      if (fc < xmax) add_crit(fc, sprintf("F = %.2f", fc), "right")
      fe <- quantile(f_vals, 0.95)
      if (fe < xmax) add_emp(fe, sprintf("sim. %.2f", fe), "right")
      
    } else if (input$dist == "chi") {
      
      chi_vals  <- chi_data()
      prev_vals <- attr(chi_vals, "prev")
      
      breaks <- pretty(range(c(chi_vals, prev_vals), na.rm = TRUE), n = 30)
      hist_with_prev(chi_vals, prev_vals, breaks,
                     main = paste("Chi-square distribution (df =", input$k_chi, ")"),
                     xlab = expression(chi^2))
      curve(dchisq(x, df = input$k_chi), add = TRUE, col = "red", lwd = 2)
      cc <- qchisq(0.95, df = input$k_chi)
      add_crit(cc, sprintf("chi-sq = %.2f", cc), "right")
      ce <- quantile(chi_vals, 0.95)
      add_emp(ce, sprintf("sim. %.2f", ce), "right")
      
    } else if (input$dist == "binom") {
      
      b_vals    <- b_data()
      prev_vals <- attr(b_vals, "prev")
      n         <- input$n_b
      p         <- input$p_b
      if (!identical(attr(prev_vals, "n"), c(n, p))) prev_vals <- NULL
      
      # One bar per possible number of successes, with the exact
      # binomial probabilities as red points
      k       <- 0:n
      p_exact <- dbinom(k, size = n, prob = p)
      breaks  <- seq(-0.5, n + 0.5, by = 1)
      hist_with_prev(b_vals, prev_vals, breaks,
                     main = paste("Binomial distribution (n =", n, ", p =", p, ")"),
                     xlab = "Number of successes",
                     ylim_extra = max(p_exact))
      lines(k, p_exact, type = "b", pch = 19, cex = 0.7, col = "red", lwd = 2)
      
      # Two-tailed rejection region (alpha = 0.05): at most 2.5% in each tail
      lo <- k[pbinom(k, n, p) <= 0.025]
      hi <- k[pbinom(k - 1, n, p, lower.tail = FALSE) <= 0.025]
      if (length(lo)) add_crit(max(lo) + 0.5, paste0("<= ", max(lo)), "left")
      if (length(hi)) add_crit(min(hi) - 0.5, paste0(">= ", min(hi)), "right")
      # Same rule applied to the simulated counts
      lo_e <- k[sapply(k, function(i) mean(b_vals <= i)) <= 0.025]
      hi_e <- k[sapply(k, function(i) mean(b_vals >= i)) <= 0.025]
      if (length(lo_e)) add_emp(max(lo_e) + 0.5, paste0("sim. <= ", max(lo_e)), "left")
      if (length(hi_e)) add_emp(min(hi_e) - 0.5, paste0("sim. >= ", min(hi_e)), "right")
      
    } else if (input$dist == "pear") {
      
      r_vals    <- r_data()
      prev_vals <- attr(r_vals, "prev")
      n         <- input$n_r
      
      # Convert each r to the t statistic that cor.test() uses
      to_t      <- function(r) r * sqrt(n - 2) / sqrt(1 - r^2)
      t_vals    <- to_t(r_vals)
      if (!is.null(prev_vals)) prev_vals <- to_t(prev_vals)
      df        <- n - 2
      
      breaks <- pretty(range(c(t_vals, prev_vals), na.rm = TRUE), n = 30)
      hist_with_prev(t_vals, prev_vals, breaks,
                     main = paste("Pearson's r converted to t (n =", n, ", df =", df, ")"),
                     xlab = expression(t == r * sqrt(n - 2) / sqrt(1 - r^2)))
      curve(dt(x, df = df), add = TRUE, col = "red", lwd = 2)
      tc <- qt(0.975, df = df)
      rc <- r_crit(n)
      add_crit(-tc, sprintf("t = %.2f (r = %.3f)", -tc, -rc), "left")
      add_crit( tc, sprintf("t = %.2f (r = %.3f)",  tc,  rc), "right")
      te <- quantile(t_vals, c(0.025, 0.975))
      add_emp(te[1], sprintf("sim. %.2f", te[1]), "left")
      add_emp(te[2], sprintf("sim. %.2f", te[2]), "right")
      
    } else if (input$dist == "spear") {
      
      sp_vals   <- sp_data()
      prev_vals <- attr(sp_vals, "prev")
      n         <- input$n_sp
      
      breaks <- seq(-1, 1, by = 0.05)
      hist_with_prev(sp_vals, prev_vals, breaks,
                     main = paste("Spearman's rank correlation when the null is true (n =", n, ")"),
                     xlab = expression(r[S]),
                     ylim_extra = dcor_null(0, n))
      curve(dcor_null(x, n), from = -0.999, to = 0.999, add = TRUE, col = "red", lwd = 2)
      rc <- r_crit(n)
      add_crit(-rc, sprintf("approx. %.3f", -rc), "left")
      add_crit( rc, sprintf("approx. %.3f",  rc), "right")
      se <- quantile(sp_vals, c(0.025, 0.975))
      add_emp(se[1], sprintf("sim. %.3f", se[1]), "left")
      add_emp(se[2], sprintf("sim. %.3f", se[2]), "right")
      
    } else {
      
      u_vals    <- u_data()
      prev_vals <- attr(u_vals, "prev")
      n1        <- input$n1_u
      n2        <- input$n2_u
      
      # U is a whole number, so use one bar per value. The axis starts well below
      # n1*n2/2 to show that U (the larger of U1 and U2) can never fall below it.
      u_low  <- floor(n1 * n2 / 4)
      breaks <- seq(u_low - 0.5, n1 * n2 + 0.5, by = 1)
      
      # The previous sample is only comparable if n1 and n2 have not changed
      if (!identical(attr(prev_vals, "n"), c(n1, n2))) prev_vals <- NULL
      
      # Exact probability of each U when the null is true. Because U is the
      # LARGER of U1 and U2, each value above the middle can arise two ways.
      u_exact <- ceiling(n1 * n2 / 2):(n1 * n2)
      p_exact <- ifelse(u_exact == n1 * n2 / 2,
                        dwilcox(u_exact, n1, n2),
                        2 * dwilcox(u_exact, n1, n2))
      
      hist_with_prev(u_vals, prev_vals, breaks,
                     main = paste("Mann-Whitney U when the null is true (n1 =", n1, ", n2 =", n2, ")"),
                     xlab = "U (the larger of U1 and U2)",
                     ylim_extra = max(p_exact))
      lines(u_exact, p_exact, type = "b", pch = 19, cex = 0.7, col = "red", lwd = 2)
      abline(v = n1 * n2 / 2, lty = 3, col = "grey50")
      text(n1 * n2 / 2, par("usr")[4] * 0.95, "n1 x n2 / 2", pos = 2, cex = 0.8, col = "grey40")
      
      # Critical value from the exact distribution (alpha = 0.05, two-tailed)
      u_crit <- min(u_exact[rev(cumsum(rev(p_exact))) <= 0.05])
      add_crit(u_crit - 0.5, paste0("U = ", u_crit), "right")
      u_all  <- u_exact
      u_sims <- u_all[sapply(u_all, function(u) mean(u_vals >= u)) <= 0.05]
      if (length(u_sims)) add_emp(min(u_sims) - 0.5, paste0("sim. U = ", min(u_sims)), "right")
    }
  })
  
  # -----------------------------
  # Explanation
  # -----------------------------
  output$explain_text <- renderUI({
    txt <- if (input$dist == "t") {
      "<b><i>t</i> distribution</b>: emerges from drawing <i>n</i> samples from a normal
       distribution and standardizing sample means:
       <i>t</i> = x&#772; &frasl; (s &frasl; &radic;<i>n</i>),
       where x&#772; is the sample mean and <i>s</i> is the sample standard deviation.
       The red line shows the theoretical <i>t</i> distribution for comparison.
       As <i>n</i> increases, the <i>t</i> distribution approaches the standard normal (Z)."
    } else if (input$dist == "f") {
      "<b><i>F</i> distribution</b>: emerges from the ratio of variances of two independent
       normal samples: <i>F</i> = s&#8321;&sup2; &frasl; s&#8322;&sup2;,
       where s&#8321;&sup2; and s&#8322;&sup2; are the variances of each group.
       Degrees of freedom: df&#8321; = n&#8321; &minus; 1 (numerator),
       df&#8322; = n&#8322; &minus; 1 (denominator).
       The red line shows the theoretical <i>F</i> distribution for comparison.
       Used extensively in ANOVA and regression to compare group variances."
    } else if (input$dist == "chi") {
      "<b>&chi;&sup2; distribution</b>: emerges from summing squared standard normal variables.
       If Z ~ N(0, 1), then &chi;&sup2;(<i>k</i>) = &sum;Z&sup2; over <i>k</i> independent draws.
       The red line shows the theoretical &chi;&sup2; distribution for comparison.
       Connects directly to variance estimation, goodness-of-fit, and contingency tables."
    } else if (input$dist == "binom") {
      "<b>Binomial distribution</b>: each simulation runs <i>n</i> independent trials, each 
       a success with probability <i>p</i> (like flipping a coin <i>n</i> times when 
       <i>p</i> = 0.5), and counts the number of successes. The blue bars show how often 
       each count occurred in the simulations; the red points are the exact probabilities 
       from the binomial formula, Pr[<i>X</i>] = <sub><i>n</i></sub>C<sub><i>X</i></sub> 
       <i>p</i><sup><i>X</i></sup>(1 &minus; <i>p</i>)<sup><i>n</i> &minus; <i>X</i></sup>, 
       calculated in R with <b>dbinom()</b>. With more simulations, the bars settle onto 
       the red points. This is the null distribution used by the binomial test."
    } else if (input$dist == "pear") {
      "<b>Pearson's correlation (<i>r</i>)</b>: each simulation draws <i>n</i> values of 
       <i>X</i> and <i>n</i> unrelated values of <i>Y</i>, so the true correlation is zero 
       (the null hypothesis is true), and calculates <i>r</i>. Pearson's correlation is 
       tested by converting <i>r</i> into a <i>t</i> statistic, 
       <i>t</i> = <i>r</i>&radic;(<i>n</i> &minus; 2) &frasl; &radic;(1 &minus; 
       <i>r</i>&sup2;), which is what <b>cor.test()</b> reports. The histogram shows these 
       converted values, and the red line is the <i>t</i> distribution with 
       <i>n</i> &minus; 2 degrees of freedom (two are lost because the regression line has 
       both a slope and an intercept). The critical values are labelled in both units: 
       as <i>t</i>, and as the equivalent value of <i>r</i>. (For the distribution of 
       <i>r</i> itself, see the Spearman option, which looks almost identical.)"
    } else if (input$dist == "spear") {
      "<b>Spearman's rank correlation (<i>r<sub>S</sub></i>)</b>: each simulation draws
       <i>n</i> values of <i>X</i> and <i>n</i> unrelated values of <i>Y</i>, so the true
       correlation is zero (the null hypothesis is true). Both variables are converted to
       ranks, and <i>r<sub>S</sub></i> is the correlation between the ranks. Even with no
       real relationship, <i>r<sub>S</sub></i> is rarely exactly zero: this histogram shows
       how far it wanders by chance. With small <i>n</i>, quite large values of
       <i>r<sub>S</sub></i> arise by chance alone, which is why the critical values in
       Whitlock &amp; Schluter's table are large for small samples. The red line is the
       null distribution of a correlation coefficient with <i>n</i> pairs (exact for
       Pearson's <i>r</i>, and a close approximation for <i>r<sub>S</sub></i>)."
    } else {
      "<b>Mann-Whitney <i>U</i></b>: each simulation draws two samples from the
       <i>same</i> normal population (the null hypothesis is true). Following Whitlock &amp;
       Schluter, all observations are ranked together, and <i>R</i><sub>1</sub> is the sum
       of the ranks for group 1. Then <i>U</i><sub>1</sub> = n&#8321;n&#8322; +
       n&#8321;(n&#8321; + 1)/2 &minus; <i>R</i><sub>1</sub>, <i>U</i><sub>2</sub> =
       n&#8321;n&#8322; &minus; <i>U</i><sub>1</sub>, and the test statistic <i>U</i> is
       the larger of the two. The red points show the exact probability of each <i>U</i>
       when the null is true. <i>U</i> can never be smaller than half of n&#8321;n&#8322;,
       and only values at or above the dashed critical value are significant. Notice that
       <i>U</i> uses only the ranks, so the same distribution applies whatever the shape
       of the original data."
    }
    crit_txt <- paste0(
      "<b>Grey dashed lines</b>: the critical value(s) at &alpha; = 0.05 (",
      if (input$dist %in% c("f", "chi")) "upper tail only" else "two-tailed",
      if (input$dist == "spear") "; for r<sub>S</sub> this is an approximation from the t distribution, and W&amp;S's table may differ slightly for small n" else "",
      "). Because every simulation here was made with the null hypothesis true, only about ",
      "5% of the simulated values should fall beyond them: those are the false positives.",
      if (input$dist %in% c("binom", "mwu"))
        " Here the values are whole numbers, so the proportion beyond the lines is at most 5%, and usually a little less." else "",
      " <span style='color:#E07B39;'><b>Orange dotted lines</b></span>: the same cut-off found ",
      "directly from the simulated values (the values that leave 5% of the simulations beyond them). ",
      "Press <b>Resample</b> and watch them bounce around the grey lines; with more simulations ",
      "they settle closer to the calculated values.")
    HTML(paste0("<p>", txt, "</p><p>", crit_txt, "</p>"))
  })
  
  output$explain_code <- renderText({
    if (input$dist == "t") {
      paste0(
        "# R code used to generate these t values\n",
        "n      <- ", input$n_t, "   # sample size\n",
        "n_sims <- ", if (!is.null(input$n_exp_t)) input$n_exp_t else 1000, "   # number of simulations\n",
        "t_vals <- replicate(n_sims, {\n",
        "  x <- rnorm(n)\n",
        "  mean(x) / (sd(x) / sqrt(n))\n",
        "})\n",
        "\n",
        "# Plot the simulated values, with the theoretical t distribution in red\n",
        "hist(t_vals, breaks = 30, probability = TRUE, col = \"#569BBD\",\n",
        "     main = \"t distribution\", xlab = \"t statistic\")\n",
        "curve(dt(x, df = n - 1), add = TRUE, col = \"red\", lwd = 2)\n",
        "abline(v = qt(c(0.025, 0.975), df = n - 1), lty = 2)   # critical values\n",
        "abline(v = quantile(t_vals, c(0.025, 0.975)), lty = 3, col = \"orange\")   # from the simulation\n"
      )
    } else if (input$dist == "f") {
      paste0(
        "# R code used to generate these F values\n",
        "n1     <- ", input$n1_f, "   # sample size group 1\n",
        "n2     <- ", input$n2_f, "   # sample size group 2\n",
        "n_sims <- ", if (!is.null(input$n_exp_f)) input$n_exp_f else 1000, "   # number of simulations\n",
        "f_vals <- replicate(n_sims, {\n",
        "  x1 <- rnorm(n1)\n",
        "  x2 <- rnorm(n2)\n",
        "  var(x1) / var(x2)\n",
        "})\n",
        "\n",
        "# Plot the simulated values (F values above 10 are not shown),\n",
        "# with the theoretical F distribution in red\n",
        "hist(f_vals, breaks = c(seq(0, 10, by = 0.25), max(f_vals, 10.25)),\n",
        "     probability = TRUE, col = \"#569BBD\", xlim = c(0, 10),\n",
        "     main = \"F distribution\", xlab = \"F statistic\")\n",
        "curve(df(x, df1 = n1 - 1, df2 = n2 - 1), add = TRUE, col = \"red\", lwd = 2)\n",
        "abline(v = qf(0.95, df1 = n1 - 1, df2 = n2 - 1), lty = 2)   # critical value\n",
        "abline(v = quantile(f_vals, 0.95), lty = 3, col = \"orange\")   # from the simulation\n"
      )
    } else if (input$dist == "chi") {
      paste0(
        "# R code used to generate these chi-square values\n",
        "k      <- ", input$k_chi, "   # degrees of freedom\n",
        "n_sims <- ", if (!is.null(input$n_exp_chi)) input$n_exp_chi else 1000, "   # number of simulations\n",
        "chi_vals <- replicate(n_sims, {\n",
        "  z <- rnorm(k)\n",
        "  sum(z^2)\n",
        "})\n",
        "\n",
        "# Plot the simulated values, with the theoretical chi-square distribution in red\n",
        "hist(chi_vals, breaks = 30, probability = TRUE, col = \"#569BBD\",\n",
        "     main = \"Chi-square distribution\", xlab = \"Chi-square\")\n",
        "curve(dchisq(x, df = k), add = TRUE, col = \"red\", lwd = 2)\n",
        "abline(v = qchisq(0.95, df = k), lty = 2)   # critical value\n",
        "abline(v = quantile(chi_vals, 0.95), lty = 3, col = \"orange\")   # from the simulation\n"
      )
    } else if (input$dist == "binom") {
      paste0(
        "# R code used to generate these binomial counts\n",
        "n      <- ", input$n_b, "   # number of trials\n",
        "p      <- ", input$p_b, "   # probability of success in each trial\n",
        "n_sims <- ", if (!is.null(input$n_exp_b)) input$n_exp_b else 1000, "   # number of simulations\n",
        "x_vals <- replicate(n_sims, sum(runif(n) < p))   # successes in n trials\n",
        "\n",
        "# Plot the simulated counts, one bar per value,\n",
        "# with the exact binomial probabilities as red points\n",
        "hist(x_vals, breaks = seq(-0.5, n + 0.5, by = 1), probability = TRUE,\n",
        "     col = \"#569BBD\", main = \"Binomial distribution\", xlab = \"Number of successes\")\n",
        "points(0:n, dbinom(0:n, size = n, prob = p), col = \"red\", pch = 19)\n"
      )
    } else if (input$dist == "pear") {
      paste0(
        "# R code used to generate these Pearson correlations\n",
        "n      <- ", input$n_r, "   # number of pairs\n",
        "n_sims <- ", if (!is.null(input$n_exp_r)) input$n_exp_r else 1000, "   # number of simulations\n",
        "r_vals <- replicate(n_sims, {\n",
        "  x <- rnorm(n)\n",
        "  y <- rnorm(n)                      # unrelated to x: null is true\n",
        "  cor(x, y)\n",
        "})\n",
        "\n",
        "# Convert each r to t, as cor.test() does\n",
        "t_vals <- r_vals * sqrt(n - 2) / sqrt(1 - r_vals^2)\n",
        "\n",
        "# Plot the converted values, with the t distribution (n - 2 df) in red\n",
        "hist(t_vals, breaks = 30, probability = TRUE, col = \"#569BBD\",\n",
        "     main = \"Pearson's r converted to t\", xlab = \"t\")\n",
        "curve(dt(x, df = n - 2), add = TRUE, col = \"red\", lwd = 2)\n",
        "abline(v = qt(c(0.025, 0.975), df = n - 2), lty = 2)   # critical values\n",
        "abline(v = quantile(t_vals, c(0.025, 0.975)), lty = 3, col = \"orange\")   # from the simulation\n",
        "\n",
        "# The critical t converted back into a critical r\n",
        "tc <- qt(0.975, df = n - 2)\n",
        "tc / sqrt(tc^2 + n - 2)\n"
      )
    } else if (input$dist == "spear") {
      paste0(
        "# R code used to generate these Spearman correlations\n",
        "n      <- ", input$n_sp, "   # number of pairs\n",
        "n_sims <- ", if (!is.null(input$n_exp_sp)) input$n_exp_sp else 1000, "   # number of simulations\n",
        "rs_vals <- replicate(n_sims, {\n",
        "  x <- rank(rnorm(n))                # rank the X values\n",
        "  y <- rank(rnorm(n))                # rank the unrelated Y values: null is true\n",
        "  cor(x, y)                          # correlation of the ranks = r_S\n",
        "})\n",
        "\n",
        "# Plot the simulated values\n",
        "hist(rs_vals, breaks = seq(-1, 1, by = 0.05), probability = TRUE,\n",
        "     col = \"#569BBD\", main = \"Spearman's rank correlation\", xlab = \"r_S\")\n",
        "\n",
        "# Proportion of simulations as extreme as a value of interest, e.g. 0.5\n",
        "mean(abs(rs_vals) >= 0.5)\n",
        "\n",
        "# Critical values of r_S from the simulation\n",
        "quantile(rs_vals, c(0.025, 0.975))\n"
      )
    } else {
      paste0(
        "# R code used to generate these Mann-Whitney U values\n",
        "n1     <- ", input$n1_u, "   # sample size group 1\n",
        "n2     <- ", input$n2_u, "   # sample size group 2\n",
        "n_sims <- ", if (!is.null(input$n_exp_u)) input$n_exp_u else 1000, "   # number of simulations\n",
        "u_vals <- replicate(n_sims, {\n",
        "  y1 <- rnorm(n1)\n",
        "  y2 <- rnorm(n2)                     # same population: null is true\n",
        "  R1 <- sum(rank(c(y1, y2))[1:n1])    # rank sum of group 1\n",
        "  U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1\n",
        "  U2 <- n1 * n2 - U1\n",
        "  max(U1, U2)                         # U = the larger of U1 and U2\n",
        "})\n",
        "\n",
        "# Plot the simulated values, one bar per value of U\n",
        "hist(u_vals, breaks = seq(min(u_vals) - 0.5, max(u_vals) + 0.5, by = 1),\n",
        "     probability = TRUE, col = \"#569BBD\",\n",
        "     main = \"Mann-Whitney U\", xlab = \"U\")\n",
        "\n",
        "# Proportion of simulations at or above a U of interest, e.g. 50\n",
        "mean(u_vals >= 50)\n"
      )
    }
  })
}

shinyApp(ui = ui, server = server)