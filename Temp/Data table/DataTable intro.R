require(data.table)
dt <- data.table(id = c(rep(1,4), rep(3,5)), a = 1:9, b=2*(1:9))

dt[, print(.SD), by = id]
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
