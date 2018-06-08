#' Run GuessMyLoi app
#'
#' @export
#'
guessApp <- function(){
  appDir <- system.file("guessMyLOIapp", package = "guessMyLOI")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `guessMyLOI`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
