source("ncv_cd.R")
source("cv_ncv_cd.R")

n <- 500
p <- 1000
p_eff <- 10
X <- matrix(rnorm(n*p, 0, 1),nrow = n)
epsilon <- rnorm(n, 0, 1)
beta <- c(runif(p_eff, 1, 2),rep(0, p-p_eff))
y <- X%*%beta + epsilon

mod_lasso <- ncv_cd(X, y, penalty = "lasso", plot = TRUE)
mod_MCP <- ncv_cd(X, y, penalty = "MCP", plot = TRUE)
mod_SCAD <- ncv_cd(X, y, penalty = "SCAD", plot = TRUE)

mod_MCP2 <- ncv_cd(X, y, penalty = "MCP", plot = TRUE, gamma = 1.2)

lasso_cv <- cv_ncv_cd(X, y, penalty = "lasso")
norm(lasso_cv$beta_cv-beta,"2")
table(pred=(lasso_cv$beta_cv!=0), true=(beta!=0))
MCP_cv <- cv_ncv_cd(X, y, penalty = "MCP")
norm(MCP_cv$beta_cv-beta,"2")
table(pred=(MCP_cv$beta_cv!=0), true=(beta!=0))
SCAD_cv <- cv_ncv_cd(X, y, penalty = "SCAD")
norm(SCAD_cv$beta_cv-beta,"2")
table(pred=(SCAD_cv$beta_cv!=0), true=(beta!=0))

mod <- cv.ncvreg(X, y, penalty = "lasso")
plot(mod)

