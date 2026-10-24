# =========================================================
# Shiny App: t vs F vs Chi-square Distributions
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
                    "Chi-square distribution" = "chi")
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
      
    } else {
      tagList(
        h4("Chi-square Distribution Settings"),
        sliderInput("k_chi", "Degrees of freedom (k):",
                    min = 1, max = 50, value = 5),
        sliderInput("n_exp_chi", "Number of simulations:",
                    min = 100, max = 5000, value = 1000, step = 100)
      )
    }
  })
  
  # -----------------------------
  # Store previous samples
  # -----------------------------
  prev_t <- reactiveVal(NULL)
  prev_f <- reactiveVal(NULL)
  prev_chi <- reactiveVal(NULL)
  
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
  
  
  # -----------------------------
  # Main plot
  # -----------------------------
  output$distPlot <- render_base("distPlot", session, function() {
    
    if (input$dist == "t") {
      
      t_vals    <- t_data()
      prev_vals <- attr(t_vals, "prev")
      df        <- input$n_t - 1
      
      breaks <- pretty(range(c(t_vals, prev_vals), na.rm = TRUE), n = 30)
      h_curr <- hist(t_vals, breaks = breaks, plot = FALSE)
      ylim   <- c(0, max(h_curr$density) * 1.15)
      
      title_txt <- paste("t distribution (n =", input$n_t, ", df =", df, ")")
      
      if (!is.null(prev_vals)) {
        hist(prev_vals, probability = TRUE, breaks = breaks,
             col = rgb(0.7, 0.7, 0.7, 0.4), border = NA,
             main = title_txt, xlab = "t statistic", ylim = ylim)
        hist(t_vals, probability = TRUE, breaks = breaks,
             col = "#569BBD", border = "black", add = TRUE)
      } else {
        hist(t_vals, probability = TRUE, breaks = breaks,
             col = "#569BBD", border = "black",
             main = title_txt, xlab = "t statistic", ylim = ylim)
      }
      
      curve(dt(x, df = df), add = TRUE, col = "red", lwd = 2)
      
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
      
    } else {
      
      chi_vals  <- chi_data()
      prev_vals <- attr(chi_vals, "prev")
      
      breaks <- pretty(range(c(chi_vals, prev_vals), na.rm = TRUE), n = 30)
      h_curr <- hist(chi_vals, breaks = breaks, plot = FALSE)
      ylim   <- c(0, max(h_curr$density) * 1.15)
      
      title_txt <- paste("Chi-square distribution (df =", input$k_chi, ")")
      
      if (!is.null(prev_vals)) {
        hist(prev_vals, probability = TRUE, breaks = breaks,
             col = rgb(0.7, 0.7, 0.7, 0.4), border = NA,
             main = title_txt, xlab = expression(chi^2), ylim = ylim)
        hist(chi_vals, probability = TRUE, breaks = breaks,
             col = "#569BBD", border = "black", add = TRUE)
      } else {
        hist(chi_vals, probability = TRUE, breaks = breaks,
             col = "#569BBD", border = "black",
             main = title_txt, xlab = expression(chi^2), ylim = ylim)
      }
      
      curve(dchisq(x, df = input$k_chi), add = TRUE, col = "red", lwd = 2)
    }
  })
  
  # -----------------------------
  # Explanation (unchanged)
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
    } else {
      "<b>&chi;&sup2; distribution</b>: emerges from summing squared standard normal variables.
       If Z ~ N(0, 1), then &chi;&sup2;(<i>k</i>) = &sum;Z&sup2; over <i>k</i> independent draws.
       The red line shows the theoretical &chi;&sup2; distribution for comparison.
       Connects directly to variance estimation, goodness-of-fit, and contingency tables."
    }
    HTML(paste0("<p>", txt, "</p>"))
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
        "})\n"
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
        "})\n"
      )
    } else {
      paste0(
        "# R code used to generate these chi-square values\n",
        "k      <- ", input$k_chi, "   # degrees of freedom\n",
        "n_sims <- ", if (!is.null(input$n_exp_chi)) input$n_exp_chi else 1000, "   # number of simulations\n",
        "chi_vals <- replicate(n_sims, {\n",
        "  z <- rnorm(k)\n",
        "  sum(z^2)\n",
        "})\n"
      )
    }
  })
}

shinyApp(ui = ui, server = server)