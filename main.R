n <- 500
p <- 1000
p_eff <- 10
X <- matrix(rnorm(n*p, 0, 1),nrow = n)
epsilon <- rnorm(n, 0, 1)
beta <- c(runif(p_eff, 1, 2),rep(0, p-p_eff))
y <- X%*%beta + epsilon
X <- scale(X)
y <- scale(y, scale = FALSE)
library(ncvreg)
mod <- ncvreg(X, y, family = "gaussian", penalty = "MCP")
plot(mod)
mod$lambda

source("ncv_cd.R")
mod <- ncv_cd(X, y, family = "gaussian", penalty = "MCP")

mod <- cv.ncvreg(X, y, family = "gaussian", penalty = "MCP")
summary(mod)

X <- scale(X)
y <- scale(y)


