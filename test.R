library(devtools)
library(ggplot2)

load_all(".")

# Check time for TA
n = 1000L
d = 10L
des = dandy::design(n, d, method = "sobol")
st = system.time({
  dandy::discrepancy(des, method = "ta", trials = 10, iters = 1000000)
})
print(st)

stop("DONE")

set.seed(1)
n = 50
methods = c("uniform", "improvedlhs", "halton", "sobol")
designs = lapply(methods, function(method) {
  dandy::design(n = n, k = 5, method = method)
})

# caluclate star discrepancy
discr.exact  = sapply(designs, dandy::discrepancy, method = "exact")
discr.approx = sapply(designs, dandy::discrepancy, method = "ta", trials = 10, iters = 1000000)

print(discr.exact)
print(discr.approx)

stop()

#set.seed(123)
x = design(n = 100, k = 4, method = "halton", u = 1, l = 0)
write.table(x, file = "~/repos/research/surrogate/projects/smboinitial/experiments/discrepancy-code/discr_calc/fuckdesign.csv", col.names = FALSE, row.names = FALSE, sep = " ")

#x = read.table("~/repos/research/surrogate/projects/smboinitial/experiments/discrepancy-code/discr_calc/fuckdesign.csv", header = TRUE)

#stop()
discr = discrepancy(x)#, method = "ta")

BBmisc::catf("Star-discr.: %.5f", discr)


stop()

methods = c("uniform", "improvedlhs", "halton", "sobol")
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
