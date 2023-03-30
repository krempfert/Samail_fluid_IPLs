#' Format number with significant digits
#' 
#' This is useful for including numbers on plots as text elements.
#' 
#' @param signif number of significant digits to display
#' @param no_sci_min use scientific notation for abs(x)<=no_sci_min, only used if always_sci = FALSE
#' @param no_sci_max use scientific notation for abs(x)>=no_sci_max, only used if always_sci = FALSE
#' @param alawys_sci set to TRUE to always use scientific notation
#' @param sci_as_latex whether to print scientific notation as latex
#' @param include_plus whether to include plus sign for positive numbers (x>0)
#' @param na_value what value to use for NA values
#' @return textual representation of the number matching the parameter settings
format_with_signif <- function(x, signif = 2, no_sci_min = 1e-5, no_sci_max = 1e5, always_sci = FALSE, sci_as_latex = FALSE, include_plus = FALSE, na_value = NA_character_) {
  if(signif <= 0) stop("this function only supports numbers with at least 1 significant digit")
  x <- base::signif(x, digits = signif)
  n_decimals = ifelse(x == 0, 0, log10(abs(x)))
  # find decimals depending on whether it's an exact power of 10 or not
  exact <- (n_decimals %% 1) < .Machine$double.eps^0.5
  n_decimals <- ifelse(exact, n_decimals, floor(n_decimals))
  # use same decimals for true 0 as the smallest abs(x)
  zeros <- abs(x) < .Machine$double.eps
  zeros[is.na(zeros)] <- FALSE
  if (any(zeros)) {
    min_decimals <- min(n_decimals[!zeros], na.rm = TRUE)
    n_decimals <- ifelse(zeros, min_decimals, n_decimals)
  }
  # determine formatting
  plus <- if(include_plus) "+" else ""
  pow <- if(sci_as_latex) "\\\\cdot{}10^{%.0f}" else "e%.0f"
  # generate output vector
  out <- character(length(x))
  out[is.na(x)] <- na_value
  # non-sci formatting
  non_sci <- !is.na(x) & !always_sci & abs(x) > no_sci_min & abs(x) < no_sci_max 
  out[non_sci] <- sprintf(sprintf("%%%s.%0.ff", plus, -n_decimals[non_sci] + signif - 1L), x[non_sci])
  # sci formatting
  sci <- !is.na(x) & !non_sci
  out[sci] <- sprintf(sprintf("%%%s.%0.ff%s", plus, signif - 1L, pow), x[sci] * 10^(-n_decimals[sci]), n_decimals[sci])
  # finished
  return(out)
}

#' export data frame to excel
#'
#' for more fine-control, call createWorkbook and add_excel_sheet directly.
#'
#' @param df data frame
#' @param df file path to the file
#' @return returns the data frame invisible for use in pipes
export_to_excel <- function(..., file, dbl_digits = 2, int_format = "0", dbl_format = sprintf(sprintf("%%.%df", dbl_digits), 0)) {
  # make excel workbook
  wb <- openxlsx::createWorkbook()
  
  # add excel shet
  sheets <- list(...)
  stopifnot(length(sheets) > 0)
  if (is.null(names(sheets)) || any(nchar(names(sheets)) == 0)) {
    names(sheets) <- paste("Sheet", seq_len(length(sheets)))
  }
  purrr::walk2(
    names(sheets),
    sheets,
    ~add_excel_sheet(wb, sheet_name = .x, df = .y, dbl_digits = dbl_digits, int_format = int_format, dbl_format = dbl_format)
  )
  
  # save workbook
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  
  # return invisible
  return(invisible(df))
}


#' add an excel sheet to a workbook
#' @param ... the data frame(s)
#' @param dbl_digits how many digits to show for dbls (all are exported)
#' @param int_format the excel formatting style for integers
#' @param dbl_format the excel formatting style for doubles (created automatically from the dbl_digits parameter)
#' @param col_max_width maximum column width
add_excel_sheet <- function(wb, sheet_name, ..., dbl_digits = 2, col_max_width = 75, int_format = "0", dbl_format = sprintf(sprintf("%%.%df", dbl_digits), 0)) {
  
  # sheet
  openxlsx::addWorksheet(wb, sheet_name)
  hs <- openxlsx::createStyle(textDecoration = "bold") # header style
  
  # data
  sheet_data_sets <- list(...)
  start_row <- 1L
  for (sheet_data in sheet_data_sets) {
    sheet_data <- dplyr::ungroup(sheet_data)
    if (ncol(sheet_data) > 0) {
      openxlsx::writeData(wb, sheet_name, sheet_data, startRow = start_row, headerStyle = hs)
      int_cols <- which(purrr::map_lgl(sheet_data, is.integer))
      dbl_cols <- setdiff(which(purrr::map_lgl(sheet_data, is.numeric)), int_cols)
      if (dbl_digits < 1) {
        int_cols <- c(int_cols, dbl_cols)
        dbl_cols <- integer()
      }
      # integer column formatting
      if (length(int_cols) > 0) {
        openxlsx::addStyle(
          wb, sheet_name, style = openxlsx::createStyle(numFmt = int_format),
          rows = (start_row + 1L):(start_row + 1L + nrow(sheet_data)),
          cols = int_cols, gridExpand = TRUE)
      }
      # double column formatting
      if (length(dbl_cols) > 0) {
        openxlsx::addStyle(
          wb, sheet_name, style = openxlsx::createStyle(numFmt = dbl_format),
          rows = (start_row + 1L):(start_row + 1L + nrow(sheet_data)),
          cols = dbl_cols, gridExpand = TRUE)
      }
      # new start row
      start_row <- start_row + nrow(sheet_data) + 2L
    }
  }
  
  # calculate header widths
  header_widths <- 
    sheet_data_sets %>% 
    # account for bold width
    purrr::map(~nchar(names(.x)))
  max_n_cols <- purrr::map_int(header_widths, length) %>% max()
  
  # calculate data widths
  if (max_n_cols > 0) {
    calculate_data_width <- function(x) {
      if (is.integer(x)) x <- sprintf("%d", x)
      else if (is.numeric(x)) x <- sprintf(paste0("%.", dbl_digits, "f"), x)
      else x <- as.character(x)
      return(max(c(0, nchar(x)), na.rm = TRUE))
    }
    data_widths <-
      sheet_data_sets %>% 
      purrr::map(
        ~dplyr::summarise_all(.x, list(calculate_data_width)) %>% 
          unlist(use.names = FALSE)
      )
    max_widths <- purrr::map2(header_widths, data_widths , ~{
      widths <- if (is.null(.y)) .x else pmax(.x, .y, 0)
      widths <- pmin(col_max_width, widths)
      c(widths, rep(0L, times = max_n_cols - length(widths)))
    })
    col_widths <- do.call(pmax, args = max_widths)
    openxlsx::setColWidths(wb, sheet_name, cols = 1:length(col_widths), widths = col_widths)
  }
  
}