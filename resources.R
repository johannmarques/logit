optim_to_latex <- function(x, label = "table:optimal_params", caption = "Optimal Parameters", param_labels = NULL) {
  require(dplyr)
  
  if (is.null(param_labels)) {
    param_labels <- paste0("Parameter ", seq_along(x$par))
  }
  
  method <- x$method
  
  table_data <- data.frame(Parameter = param_labels,
                           Value = round(x$par,3))
  
  latex_code <- paste0("\\begin{table}[htb]\n",
                       "\\centering\n",
                       "\\caption{", caption, "}\\label{", label, "}\n",
                       "\\begin{tabular}{lc}\n \\hline \n")
  
  for (i in seq(nrow(table_data))) {
    latex_code <- paste0(latex_code, table_data[i, "Parameter"], " & ", table_data[i, "Value"], " \\\\ \n")
  }
  
  latex_code <- paste0(latex_code,
                       "\\hline\n", "\\end{tabular}\n",
                       "\\end{table}\n")
  
  cat(latex_code)
  return(latex_code)
}

my_col <- c(
  '#274C77',
           '#6096BA',
           '#A3CEF1',
           '#8B8C89'
)

my_palette <- function(n = 1){
  return(my_col[1:n])
}

my_theme <- theme_light(base_size = 20) +
  theme(
    strip.background=element_rect(colour=my_palette(1),
                                  fill=my_palette(1)),
    axis.text = element_text(color = 'black'),
    text = element_text(family = 'Lato')
  )
