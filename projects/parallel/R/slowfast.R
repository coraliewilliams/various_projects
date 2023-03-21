####### Make my code faster
# Take variables out of the dataframe and update them one by one


###### SLOW CODE - 102 seconds

n = 60000
r = 10000
d = data.frame(w = numeric(n), x = numeric(n),
               y = numeric(n), z = numeric(n))
value = c(1, 2, 3, 4)
system.time({
  for(i in 1:r) {
    j = sample(n, 1)
    d[j,] = value
  }
})

######## FAST CODE - 0.23 seconds

n = 60000
r = 10000
w = x = y = z = numeric(n)
value = c(1, 2, 3, 4)
system.time({
  for(i in 1:r) {
    j = sample(n, 1)
    w[j] = value[1]
    x[j] = value[2]
    y[j] = value[3]
    z[j] = value[4]
  }
})



