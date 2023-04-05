set.seed(12)

# Definimos una matriz aleatoria A de tamaño 10x5
A <- matrix(rnorm(50), ncol = 5)

# Definimos un vector aleatorio b de longitud 10
set.seed(11)
b <- rnorm(10)

# Definimos un vector x de longitud 5
set.seed(13)
x <- rnorm(5)

# Definimos un valor de eta
eta <- 0.1

# Ejecutamos la función fw.step
fw.step(A, x, b, eta)

alpha = NULL



Ax = A %*% x

half.grad = t(Ax - b) %*% A + eta * x
i = which.min(half.grad)
if (!is.null(alpha)) {
  x = x * (1 - alpha)
  x[i] = x[i] + alpha
  return(x)
} else {
  d.x = -x; d.x[i] = 1 - x[i]
  if (all(d.x == 0)) { return(x) }
  d.err = A[, i] - Ax
  step = -t(c(half.grad)) %*% d.x / (sum(d.err^2) + eta * sum(d.x^2))
  constrained.step = min(1, max(0, step))
  return(x + constrained.step * d.x)
}

