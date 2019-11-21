#
# This is a Shiny web application to play around with the logic behind Support Vector Machines.
#
#

library(plotly)
library(purrr)
library(shiny)
library(matlib)
library(SVMMaj)

# Define UI for the application
ui <- fluidPage(
  titlePanel("SVM intuition plotter"),
  
  fluidRow(
    column(8, plotlyOutput("p")),
    column(4, textOutput("curr_loss"))
    ),
    br(),
    br(),
  fluidRow(
    column(2, numericInput("marg", "Margin:", 0.1, min = 0.1, max = 10, step = 0.1)),
    column(2, numericInput("intercept", "Intercept:", 0, min = -10, max = 10, step = 0.1)),
    column(2, numericInput("lambda", "Lambda:", 1, min = 0, max = 1000)),
    column(3, actionButton("showproj", "Show/Hide projected points")),
    column(3, actionButton("showsol", "Show/Hide optimal solution")),
    #,column(2, numericInput("angle", "Angle:", 1, min = 0, max = 1000))
    br()
    ),
  fluidRow(
    column(2, actionButton("genpoints", "Generate new sample"))
    ),
    br(),
  fluidRow(  
    column(6, HTML(
      paste(
        h3("Instructions:"),#'<br/>',
        h4(tags$b("KNOWN BUG:"), "To make sure everything updates correctly, move the blue dot first."),
        h4("Use the mouse to move the blue dot to change the orientation of the axes."),
        h4("You can set the lambda, margin and intercept using the input boxes. You can toggle the 
           optimal solution for this lambda as found by", tags$a(href = 'https://cran.rstudio.com/web/packages/SVMMaj/index.html', 'SVMMaj')," package."),
        h4(tags$b("Disclaimer:"), " This was made as a quick proof of concept/experiment and to explore Shiny and plotly. Suggestions and bug reports welcome! ")
      )
     )
    ),
    column(6, verbatimTextOutput("curr_loss2"))
    )
)
# Calculate loss function
loss <- function(xy, type, wx, wy, c_intercept, lambda=1){
  sum(pmax(0, (abs(type) - sign(type) * (c_intercept + xy$x * wx + xy$y * wy )))) + lambda * (wx^2 + wy^2)
}

# Generate a new sample of points
generate_points <- function(n=20, meanx=c(-2, 2), meany=c(-1, 1)){
  list(x = c(rnorm(n, mean = meanx[1]), rnorm(n, mean = meanx[2])),
       y = c(rnorm(n, mean = meany[1]), rnorm(20, mean = meany[2])),
       type = c(rep(-1, n), rep(1, n)),
       type_color = c(rep(I('red'), n), rep(I('green'), n)))
}

# Project points on axes w1 and w2
proj_points <- function(x, y, w1, w2){
  print(w1)
  print(w2)
  map2_dfr(c(x = x), c(y = y), function(x, y, ...){
    output <- Proj(c(x, y), ...);
    list(x = output[1], y = output[2])}, 
    X = c(w1, w2))
}

# Calculate new coords for axes 
set_axes <- function(new_pos, w, marg, intercept){
  ang <- angle(c(0, 1), c(new_pos), degree = FALSE)
  ang <- ifelse(new_pos[1] > 0 , ang, -ang)
  w$y1 <- c( 2 * sin(ang), -2 * sin(ang)) + intercept * cos(ang)
  w$x1 <- c(-2 * cos(ang),  2 * cos(ang)) + intercept * sin(ang)
  
  w$y2 <- c(-2 * cos(ang), 2 * cos(ang))
  w$x2 <- c(-2 * sin(ang), 2 * sin(ang))
  
  w$ym1 <- w$y1 + (marg) * cos(ang)
  w$ym2 <- w$y1 - (marg) * cos(ang)
  
  w$xm1 <- w$x1 + (marg) * sin(ang)
  w$xm2 <- w$x1 - (marg) * sin(ang)
  w
}


server <- function(input, output, session) {
  # It is necessary to initialise these values to setup projected points
  points <- generate_points()
  init_x <- points$x
  init_y <- points$y
  type <- points$type
  type_color <- points$type_color
  
  # The original points
  xy <- reactiveValues(
    x = init_x,
    y = init_y
  )
  
  # Projected points  
  proj_xy <- reactiveValues(
    x = rep(0, length(init_x)),
    y = init_y,
    color = type_color,
    opacity = 0
  )
  
  cursor_coord <- reactiveValues(
    x = 0,
    y = 3 
  )
  
  # Axes + margins for plotting
  waxes <- reactiveValues(
    x1 = c(-2, 2),
    y1 = c(0, 0),
    x2 = c(0, 0),
    y2 =  c(-2, 2),
    xm1 = c(-2, 2),
    ym1 = c(-0.1, -0.1),
    xm2 = c(-2, 2),
    ym2 = c(0.1, 0.1)
  )
  
  # Optimal solution for this lambda
  solution <- svmmaj(cbind(init_x, init_y), type, lambda = 1, scale = "none")
  # x and y for solution (not normalised)
  new_pos <- SVMMaj:::beta.theta(solution$method, solution$theta)
  # margin for solution
  marg_opt <- sum(new_pos^2)/2
  
  #initialise temp variable for use in reactiveValues
  wtemp <- set_axes(new_pos, list(x1 = 0, x2 = 0, y1 = 0, y2 = 0, 
                                  xm1 = 0, xm2 = 0, ym1 = 0, ym2 = 0),
                    marg_opt, solution$theta[1])
  
  # Axes + margins for the solution
  wsols <- reactiveValues(
    x1 = wtemp$x1, 
    y1 = wtemp$y1, 
    x2 = wtemp$x2, 
    y2 = wtemp$y2, 
    xm1 = wtemp$xm1, 
    ym1 = wtemp$ym1, 
    xm2 = wtemp$xm2, 
    ym2 = wtemp$ym2, 
    opacity = 0
  )

  # Render graph
  output$p <- renderPlotly({
    # shape for the cursor
    cursor <- list(
      type = "circle",
      xanchor = cursor_coord$x,
      yanchor = cursor_coord$y,
      # diameter is 2 pixels
      x0 = -4, x1 = 4,
      y0 = -4, y1 = 4,
      xsizemode = "pixel", 
      ysizemode = "pixel",
      fillcolor = "blue",
      line = list(color = "transparent")
    )
    
    
    # plot the poinst and lines
    p <- plot_ly(type = 'scatter', mode = 'markers') %>%
      # plot axes
      add_lines(x = waxes$x1, y = waxes$y1, color = I("green"), mode = "lines") %>%
      add_lines(x = waxes$x2, y = waxes$y2, color = I("black"), mode = "lines") %>%
      # Plot margins
      add_lines(x = waxes$xm1, y = waxes$ym1, color = I("gray"), mode = "lines") %>%
      add_lines(x = waxes$xm2, y = waxes$ym2, color = I("gray"), mode = "lines") %>%
      
      # Plot axes and margins of the optimal solution
      #ldash <- list(dash = 'dash') # Gives error if used in lines below
      add_lines(x = wsols$x1, y = wsols$y1, color = I("magenta"), opacity = wsols$opacity, mode = "lines", line = list(dash = 'dash')) %>%
      add_lines(x = wsols$x2, y = wsols$y2, color = I("magenta"), opacity = wsols$opacity, mode = "lines", line = list(dash = 'dash')) %>%
      add_lines(x = wsols$xm1, y = wsols$ym1, color = I("magenta"), opacity = wsols$opacity, mode = "lines", line = list(dash = 'dash')) %>%
      add_lines(x = wsols$xm2, y = wsols$ym2, color = I("magenta"), opacity = wsols$opacity, mode = "lines", line = list(dash = 'dash')) %>%
      
      # Plot the original points
      add_trace(x = xy$x, y = xy$y, mode = "markers", marker = list(color = type_color)) %>% 
      # Plot the points projected on the axes
      add_trace(x = proj_xy$x, y = proj_xy$y, mode = 'markers', 
                marker = list(
                  color = 'rgba(135, 206, 250, 0)', # Set fill transparent
                  opacity = proj_xy$opacity,
                  size = 5,
                  line = list(
                    color = type_color,
                    width = 1
                  )
                )
              ) %>% 
      # Plot cursor as shape
      layout(shapes = cursor, 
             xaxis = list(
               fixedrange = TRUE, # Needed to prevent zooming
               #range = c(-4, 4),
               scaleanchor = 'y', # Equal aspect ratio
               scaleratio = 1, 
               constrain = "domain" ),
             yaxis = list(
               range = c(-4, 4),
               fixedrange = TRUE
               ),
             showlegend = FALSE,
             autosize = F 
             ) %>%
      config(edits = list(shapePosition = TRUE), # No way to disable resize on shapes
             displayModeBar = FALSE, 
             scrollZoom = FALSE 
             )
  })
  
  output$curr_loss <- renderPrint({
    cat("Current loss: ", 0)
  })
  
  output$curr_loss2 <- renderPrint({
    cat("Click Show solution to show the parameters of \n the optimal solution in this window")
  })
  
  output$opt_loss <- renderPrint({
    cat("Optimal (solution) loss: ", 0)
  })
  
  # Switch opacity of solution
  observeEvent(input$showsol, {
    wsols$opacity = (wsols$opacity + 1) %% 2
  })
  
  # Switch opacity of projected points
  observeEvent(input$showproj, {
    proj_xy$opacity = (proj_xy$opacity + 1) %% 2
  })
  
  # Generate new points
  observeEvent(input$genpoints, {
    points <- generate_points()
    xy$x <- points$x
    xy$y <- points$y
  })
  
  # update all reactive values when cursor is moved
  observe({
    ed <- event_data("plotly_relayout")
    shape_anchors <- ed[grepl("^shapes.*anchor$", names(ed))]
    if (length(shape_anchors) != 2) return()

    new_pos <- c(shape_anchors[[1]], shape_anchors[[2]])
    # Calculate angle from y axis 
    ang <- angle(c(0, 1), c(new_pos), degree = FALSE)
    ang <- ifelse(new_pos[1] > 0 , ang, -ang)
    
    # Sanity check marg
    curr_marg <- input$marg
    if (is.na(curr_marg) | curr_marg == 0) { 
      curr_marg <- 0.0001
    }
    
    new_pos <- curr_marg * new_pos
    cursor_coord$x <- 1.75 * new_pos[1]
    cursor_coord$y <- 1.75 * new_pos[2]

    set_axes(new_pos, waxes, curr_marg, input$intercept)
    
    # Calculate and set projected points
    proj_ax <- proj_points(xy$x, xy$y, waxes$x2[2], waxes$y2[2])
    proj_xy$y <- proj_ax$y
    proj_xy$x <- proj_ax$x
    
    # Sanity check lambda
    curr_lambda <- input$lambda
    if (is.na(curr_lambda)) { 
      curr_lambda <- 0
    }
    
    closs <- loss(xy, type, new_pos[1], new_pos[2], input$intercept, lambda = curr_lambda)
    output$curr_loss <- renderText(paste("Current loss:", round(closs, digits = 4)))
    if (wsols$opacity > 0) {
      output$curr_loss2 <- renderText(paste("Optimal model:\n",
                                            "\nLoss:",round(solution$loss, digits = 4),
                                            "\nx,y:", round(new_pos[1], digits = 4), round(new_pos[2], digits = 4),
                                            "\nmargins:", round(marg_opt, digits = 4),
                                            "\nintercept:", round(solution$theta[1], digits = 4) )
                                      )
    } else {
      output$curr_loss2 <- renderText(paste("Click Show solution to show the parameters of \n the optimal solution in this window"))
    }
    
    solution <- svmmaj(cbind(xy$x, xy$y), type, lambda = curr_lambda, scale = "none")
    new_pos <- SVMMaj:::beta.theta(solution$method, solution$theta)
    marg_opt <- sum(new_pos^2)/2
    
    # Change solution axes
    wsols <- set_axes(new_pos, wsols, marg_opt, solution$theta[1])

  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

