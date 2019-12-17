library(devtools)
library(lhs)
library(ggplot2)


load_all(".")

methods = c("uniform", "lhs", "halton", "sobol")
reps = 1:2L
gr = expand.grid(methods, reps)

designs = apply(gr, 1L, function(exp) {
  method = exp[1L]; repl = exp[2L]
  des = design(n = 100, k = 2, l = -5, u = c(5, 10), method = method)
  des$method = method
  des$repl = repl
  return(des)
})
designs = do.call(rbind, designs)

pl = ggplot(data = designs, aes(x = x1, y = x2))
pl = pl + geom_point()
pl = pl + facet_grid(repl ~ method, labeller = ggplot2::label_both)
pl = pl + theme_bw()
print(pl)

x = design(n = 100, k = 2, method = "halton")

#write.table(x, file = "lhs.csv", col.names = FALSE, row.names = FALSE, sep = " ")


discr = stardiscrepancy(x)

BBmisc::catf("Star-discr.: %.5f", discr)
