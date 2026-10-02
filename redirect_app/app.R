# Lightweight redirect app for the legacy shinyapps.io URL.
# Replaces the original HOHC app once the new site is confirmed stable.
# The published URL in the paper keeps working: visitors are forwarded
# to the new permanent home of HOHC.

library(shiny)

NEW_URL <- "https://hohc.pages.dev/"

ui <- fluidPage(
  tags$head(
    tags$title("HOHC has moved"),
    tags$script(HTML(sprintf('window.location.replace("%s");', NEW_URL)))
  ),
  div(
    style = "font-family: sans-serif; max-width: 600px; margin: 80px auto; text-align: center;",
    h2("HOHC has moved"),
    p("The Human Organ Hormonal Communication atlas is now hosted at:"),
    p(a(href = NEW_URL, NEW_URL, style = "font-size: 1.2em;")),
    p("You will be redirected automatically.")
  )
)

server <- function(input, output, session) {}

shinyApp(ui, server)
