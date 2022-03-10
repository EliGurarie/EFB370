# lifted entirely from Rage::plot_life_cycle
leslie_diagram <- function (M, stages, title = NULL, shape = "egg", fontsize = 10, 
          nodefontsize = 12, edgecol = "grey", ...) 
{
  if (missing(stages) && is.null(dimnames(M))) {
    stages <- seq_len(ncol(M))
  }
  else if (missing(stages) && !is.null(dimnames(M))) {
    stages <- dimnames(M)[[1]]
    if (!identical(dimnames(M)[[1]], dimnames(M)[[2]])) {
      message(strwrap(prefix = " ", initial = "", "Dimension names of 'M' are not identical \n                for rows and columns. Using row names."))
    }
  }
  
  graph <- expand.grid(to = stages, from = stages)
  graph$trans <- round(c(M), 3)
  graph <- graph[graph$trans > 0, ]
  nodes <- paste(paste0("'", stages, "'"), collapse = "; ")
  graph$min_len <- (as.numeric(graph$to) - as.numeric(graph$from)) * 3
  edges <- paste0("'", graph$from, "'", " -> ", "'", graph$to, 
                  "'", "[minlen=", graph$min_len, ",fontsize=", fontsize, 
                  ",color=", edgecol, ",xlabel=", paste("\"", graph$trans), 
                  "\"]\n", collapse = "")
  grViz(paste("\ndigraph {\n  {\n    graph[overlap=false];\n    rank=same;\n    node [shape=", 
              shape, ", fontsize=", nodefontsize, "];", nodes, " \n  }", 
              "ordering=out\n  x [style=invis]\n  x -> {", nodes, "} [style=invis]", 
              edges, "labelloc=\"t\";\n  label=\"", title, "\"\n}"), ...)
}