library(shiny)
library(ggplot2)

# Define UI
ui <- fluidPage(
  titlePanel("Numerical Integration Methods"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("func", "Function f(x):", "sin(x)"),
      numericInput("a", "Lower limit (a):", 0),
      numericInput("b", "Upper limit (b):", pi),
      numericInput("n", "Number of intervals (n):", 10, min = 1),
      selectInput("method", "Integration Method:",
                  choices = c("Riemann", "Midpoint", "Trapezoid", "Simpson"))
    ),
    
    mainPanel(
      verbatimTextOutput("values"),
      plotOutput("plot", height = "500px"),
      uiOutput("formula"),
      uiOutput("description")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive function expression
  f <- reactive({
    function(x) eval(parse(text = input$func), list(x = x))
  })
  
  # Compute integration values
  integration_values <- reactive({
    a <- input$a
    b <- input$b
    n <- input$n
    h <- (b - a) / n
    x <- seq(a, b, length.out = n + 1)
    mid_x <- (x[-1] + x[-length(x)]) / 2
    fx <- f()
    
    exact <- integrate(f(), a, b)$value
    estimate <- switch(input$method,
                       "Riemann" = sum(fx(x[-length(x)]) * h),
                       "Midpoint" = sum(fx(mid_x) * h),
                       "Trapezoid" = sum((fx(x[-1]) + fx(x[-length(x)])) * h / 2),
                       "Simpson" = {
                         if (n %% 2 == 1) n <- n + 1  # Simpson requires even n
                         h <- (b - a) / n
                         x_sim <- seq(a, b, length.out = n + 1)
                         fx_sim <- fx(x_sim)
                         (h / 3) * (fx_sim[1] + 2 * sum(fx_sim[seq(3, n-1, by = 2)]) +
                                      4 * sum(fx_sim[seq(2, n, by = 2)]) + fx_sim[n + 1])
                       }
    )
    
    list(exact = exact, estimate = estimate, x = x, fx = fx)
  })
  
  # Output integration values
  output$values <- renderPrint({
    vals <- integration_values()
    cat("Exact Integral:", vals$exact, "\n")
    cat("Estimated Value:", vals$estimate)
  })
  
  # Plot the approximation
  output$plot <- renderPlot({
    vals <- integration_values()
    a <- input$a
    b <- input$b
    n <- input$n
    h <- (b - a) / n
    x <- vals$x
    fx <- vals$fx
    
    x_range <- seq(a, b, length.out = 1000)
    df_func <- data.frame(x = x_range, y = fx(x_range))
    
    base_plot <- ggplot() +
      geom_line(data = df_func, aes(x = x, y = y), color = "red", size = 1.2) +
      labs(title = paste(input$method, "Approximation"), x = "x", y = "f(x)") +
      xlim(a, b)
    
    if (input$method == "Riemann") {
      df <- data.frame(
        xmin = x[-length(x)],
        xmax = x[-1],
        ymin = 0,
        ymax = fx(x[-length(x)])
      )
      base_plot <- base_plot +
        geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "lightblue", color = "blue", alpha = 0.5, inherit.aes = FALSE)
      
    } else if (input$method == "Midpoint") {
      mid_x <- (x[-1] + x[-length(x)]) / 2
      df <- data.frame(
        xmin = x[-length(x)],
        xmax = x[-1],
        ymin = 0,
        ymax = fx(mid_x)
      )
      base_plot <- base_plot +
        geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "orange", color = "darkorange", alpha = 0.5, inherit.aes = FALSE)
      
    } else if (input$method == "Trapezoid") {
      trap_data <- do.call(rbind, lapply(1:n, function(i) {
        data.frame(
          x = c(x[i], x[i+1], x[i+1], x[i]),
          y = c(0, 0, fx(x[i+1]), fx(x[i])),
          group = i
        )
      }))
      base_plot <- base_plot +
        geom_polygon(data = trap_data, aes(x = x, y = y, group = group),
                     fill = "lightgreen", color = "darkgreen", alpha = 0.4, inherit.aes = FALSE)
      
    } else if (input$method == "Simpson") {
      if (n %% 2 == 1) {
        n <- n + 1
        x <- seq(a, b, length.out = n + 1)
      }
      
      for (i in seq(1, n, by = 2)) {
        x_left <- x[i]
        x_mid <- (x[i] + x[i + 1]) / 2
        x_right <- x[i + 1]
        
        y_left <- fx(x_left)
        y_mid <- fx(x_mid)
        y_right <- fx(x_right)
        
        mat <- matrix(c(x_left^2, x_left, 1,
                        x_mid^2, x_mid, 1,
                        x_right^2, x_right, 1), nrow = 3, byrow = TRUE)
        y_vals <- c(y_left, y_mid, y_right)
        
        coeffs <- solve(mat, y_vals)
        a_q <- coeffs[1]; b_q <- coeffs[2]; c_q <- coeffs[3]
        
        x_parabola <- seq(x_left, x_right, length.out = 100)
        y_parabola <- a_q * x_parabola^2 + b_q * x_parabola + c_q
        
        parabola_df <- data.frame(x = x_parabola, y = y_parabola)
        
        base_plot <- base_plot +
          geom_line(data = parabola_df, aes(x = x, y = y), color = "blue", size = 1.2) +
          geom_area(data = parabola_df, aes(x = x, y = y), fill = "purple", alpha = 0.4)
      }
    }
    
    base_plot
  })
  
  # Latex formula
  output$formula <- renderUI({
    eq <- switch(input$method,
                 "Riemann" = "$$\\int_a^b f(x)\\,dx \\approx \\sum_{i=0}^{n-1} f(x_i) \\cdot h$$",
                 "Midpoint" = "$$\\int_a^b f(x)\\,dx \\approx \\sum_{i=0}^{n-1} f\\left(\\frac{x_i + x_{i+1}}{2}\\right) \\cdot h$$",
                 "Trapezoid" = "$$\\int_a^b f(x)\\,dx \\approx \\sum_{i=0}^{n-1} \\frac{f(x_i) + f(x_{i+1})}{2} \\cdot h$$",
                 "Simpson" = "$$\\int_a^b f(x)\\,dx \\approx \\frac{h}{3}\\left[f(x_0) + 4 \\sum_{\\text{odd } i=1}^{n-1} f(x_i) + 2 \\sum_{\\text{even } i=2}^{n-2} f(x_i) + f(x_n)\\right]$$"
    )
    withMathJax(HTML(paste("<h4>Formula:</h4>", eq)))
  })
  
  # Explanation
  output$description <- renderUI({
    desc <- switch(input$method,
                   "Riemann" = "The Riemann sum approximates the integral by using the left endpoints of subintervals to construct rectangles under the curve.",
                   "Midpoint" = "The Midpoint rule uses the function value at the middle of each interval to estimate the area under the curve.",
                   "Trapezoid" = "The Trapezoid rule approximates the region under the curve as a series of trapezoids, offering better accuracy than rectangles.",
                   "Simpson" = "Simpson's rule uses parabolas to approximate the curve over pairs of subintervals. It requires an even number of intervals and gives higher accuracy."
    )
    HTML(paste("<h4>Explanation:</h4><p>", desc, "</p>"))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
