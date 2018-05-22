require(data.table)
dt <- data.table(id = c(rep(1,4), rep(3,5)), a = 1:9, b=2*(1:9))
setkey(dt, id)

dt[, TstarById(.SD), by = id] # <-- this

idvec <- c(1,3)

for(i in idvec){
  TstarById(dt[dt[["id"]]==i,]) # <-- equivalent
  
}

dt[, print(.SD), by = id, .SDcols = "a"]
dt[, print(.SD[["a"]]), by = id]
dt[, print(.N), by = id]
dt[, print(.BY), by = id]
dt[, print(.SD[.N,]), by = id]

# Using a function that takes a vector
dt[, sum(.SD), by = id, .SDcols = "a"] #Error read msg!  Issue: Sum cannot be applied to a data.table
dt[, sum(.SD[["a"]]), by = id]         # OK
dt[, lapply(.SD, sum), by = id, .SDcols = c("a")] # OK

dt[, lapply(.SD, sum), by = id, .SDcols = c("a", "b")] #Another OK
dt[, lapply(.SD, sum), by = id] # Same as above

testfun <- function(x,y) return(x*y)
dt[, testfun(a,b), by = .I]
