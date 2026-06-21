# =========================================================
# Shiny App: t vs F vs Chi-square Distributions
# (with controlled resampling + previous sample overlay)
# =========================================================

# https://hbctraining.github.io/Training-modules/RShiny/lessons/shinylive.html
# Run the shinylive::export line to populate the docs folder 
# so that shinylive works from github
#shinylive::export(appdir = "../DerivingClassicDistributions/", destdir = "docs")
#httpuv::runStaticServer("docs/", port = 8008)

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

      br(),
    
      ),
    
   
    
    mainPanel(
      plotOutput("distPlot", height = "500px"),
      verbatimTextOutput("explain")
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
  output$distPlot <- renderPlot({
    
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
      
      breaks <- pretty(range(c(f_vals, prev_vals), na.rm = TRUE), n = 30)
      h_curr <- hist(f_vals, breaks = breaks, plot = FALSE)
      ylim   <- c(0, max(h_curr$density) * 1.15)
      
      title_txt <- paste("F distribution (df1 =", df1, ", df2 =", df2, ")")
      
      if (!is.null(prev_vals)) {
        hist(prev_vals, probability = TRUE, breaks = breaks,
             col = rgb(0.7, 0.7, 0.7, 0.4), border = NA,
             main = title_txt, xlab = "F statistic",
             xlim = c(0, 10), ylim = ylim)
        hist(f_vals, probability = TRUE, breaks = breaks,
             col = "#569BBD", border = "black", add = TRUE, xlim = c(0, 10))
      } else {
        hist(f_vals, probability = TRUE, breaks = breaks,
             col = "#569BBD", border = "black",
             main = title_txt, xlab = "F statistic",
             xlim = c(0, 10), ylim = ylim)
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
  output$explain <- renderText({
    
    if (input$dist == "t") {
      paste0(
        "t distribution: emerges from drawing n samples from a normal distribution\n",
        "and standardizing sample means: t = mean(x) / (sd(x) / sqrt(n)).\n",
        "Red line shows the theoretical t distribution for comparison.\n",
        "As n increases, the t distribution approaches the standard normal.\n\n",
        "--- R code used to generate these t values ---\n",
        "n      <- ", input$n_t, "   # sample size\n",
        "n_sims <- ", if (!is.null(input$n_exp_t)) input$n_exp_t else 1000, "   # number of simulations\n",
        "t_vals <- replicate(n_sims, {\n",
        "  x <- rnorm(n)\n",
        "  mean(x) / (sd(x) / sqrt(n))\n",
        "})\n"
      )
      
    } else if (input$dist == "f") {
      paste0(
        "F distribution: emerges from the ratio of variances of two independent normal samples.\n",
        "Degrees of freedom: df1 = n1-1 (numerator), df2 = n2-1 (denominator).\n",
        "Red line shows the theoretical F distribution for comparison.\n",
        "Used extensively in ANOVA and regression to compare group variances.\n\n",
        "--- R code used to generate these F values ---\n",
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
        "Chi-square distribution: emerges from summing squared standard normal variables.\n",
        "If Z ~ N(0,1), then chi-sq(k) = sum(Z^2) over k independent draws.\n",
        "Red line shows the theoretical chi-square distribution for comparison.\n",
        "Connects directly to variance estimation, goodness-of-fit, and contingency tables.\n\n",
        "--- R code used to generate these chi-square values ---\n",
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