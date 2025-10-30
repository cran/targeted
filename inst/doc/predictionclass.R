## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
 collapse = TRUE,
 #dev="png",
 comment = "#>"
)
library("targeted")

## ----pbcdata------------------------------------------------------------------
data(pbc, package="survival")
pbc <- transform(pbc, y = (time < 730) * (status > 0))

## ----logistic1----------------------------------------------------------------
lr <- learner_glm(y ~ age, family = binomial())
lr$estimate(pbc)

## ----logistic2----------------------------------------------------------------
lr$predict(newdata = data.frame(age = c(20, 40, 60, 80)))

## ----eval = FALSE-------------------------------------------------------------
# ?learner # help(learner)

## -----------------------------------------------------------------------------
lr_sl <- learner_sl(
  learners = list(
    glm1 = learner_glm(y ~ age * bili, family = "binomial"),
    glm2 = learner_glm(y ~ age, family = "binomial"),
    gam = learner_gam(y ~ s(age) + s(bili), family = "binomial")
  )
)
lr_sl$estimate(pbc, nfolds = 10)
lr_sl

## ----learner_expand_grid------------------------------------------------------
lrs <- learner_expand_grid( learner_xgboost,
                            list(formula = Sepal.Length ~ .,
                                eta = c(0.2, 0.5, 0.3)) )
lrs

## -----------------------------------------------------------------------------
lr_xgboost <- learner_xgboost(
  formula = y ~ age + sex + bili,
  eta = 0.3, nrounds = 5,  # hyperparameters
  objective = "binary:logistic" # learning task
)
lr_xgboost

## -----------------------------------------------------------------------------
lr_xgboost$estimate(data = pbc)

## -----------------------------------------------------------------------------
class(lr_xgboost$fit)

## -----------------------------------------------------------------------------
lr_xgboost$predict(head(pbc))

## -----------------------------------------------------------------------------
lr <- learner_glm(y ~ age, family = "binomial")
estimate(lr, pbc)
predict(lr, head(pbc))

## -----------------------------------------------------------------------------
lr_xgboost$summary()$estimate

## -----------------------------------------------------------------------------
lr_xgboost$update(y ~ age + sex)

## -----------------------------------------------------------------------------
head(lr_xgboost$design(pbc)$x)

## -----------------------------------------------------------------------------
head(lr_xgboost$response(pbc))

## -----------------------------------------------------------------------------
# future::plan("multicore")
lrs <- list(
  glm = learner_glm(y ~ age + age, family = "binomial"),
  gam = learner_gam(y ~ s(age) + s(bili), family = "binomial")
)
 # 2 times repeated 5-fold cross-validation
cv(lrs, data = pbc, rep = 2, nfolds = 5)

## -----------------------------------------------------------------------------
new_glm <- function(formula, ...) {
  learner$new(
    formula = formula,
    estimate = stats::glm,
    predict = stats::predict,
    predict.args = list(type = "response"),
    estimate.args = list(...),
    info = "new glm learner" # optional
  )
}
lr <- new_glm(y ~ age, family = "binomial")
lr

## -----------------------------------------------------------------------------
lr$estimate(pbc)

## -----------------------------------------------------------------------------
fit <- glm(y ~ age, family = "binomial", data = pbc)
all(coef(fit) == coef(lr$fit))

## -----------------------------------------------------------------------------
lr$predict(head(pbc))

## -----------------------------------------------------------------------------
predict(fit, newdata = head(pbc), type = "response")

## -----------------------------------------------------------------------------
lr$estimate(pbc, family = "poisson")
lr$predict(head(pbc), type = "link")

## -----------------------------------------------------------------------------
new_grf <- function(formula, ...) {
  learner$new(
    formula = formula,
    estimate = function(x, y, ...) grf::probability_forest(X = x, Y = y, ...),
    predict = function(object, newdata) {
      predict(object, newdata = newdata)$predictions
    },
    estimate.args = list(...),
    info = "grf::probability_forest"
  )
}
lr <- new_grf(as.factor(y) ~ age + bili, num.trees = 100)
lr$estimate(pbc)

## -----------------------------------------------------------------------------
dsgn <- lr$design(pbc)

## -----------------------------------------------------------------------------
new_sl <- function(learners, ...) {
  learner$new(
    info = "new superlearner",
    estimate = superlearner,
    predict = targeted:::predict.superlearner,
    estimate.args = c(list(learners = learners), list(...))
  )
}
lrs <- list(
  glm = learner_glm(y ~ age, family = "binomial"),
  gam = learner_gam(y ~ s(age), family = "binomial")
)
lr <- new_sl(lrs, nfolds = 2)
lr$estimate(pbc)
lr

## ----naivebayes_aggregate_data------------------------------------------------
library("data.table")
dd <- data.table(pbc)[!is.na(trt), .(.N), by=.(y,trt,sex)]
print(dd)

## ----naivevbayes--------------------------------------------------------------
lr <- learner_naivebayes(y ~ trt + sex + weights(N))
lr$estimate(dd)
lr$predict(dd)

## -----------------------------------------------------------------------------
targeted:::weights.numeric

## ----eststrata----------------------------------------------------------------
est <- function(formula, data, strata, ...)
  lapply(levels(strata), \(s) lm(formula, data[which(strata==s),]))

pred <- function(object, newdata, strata, ...) {
  res <- numeric(length(strata))
  for (i in seq_along(levels(strata))) {
    idx <- which(strata == levels(strata)[i])
    res[idx] <- predict(object[[i]], newdata[idx, ], ...)
  }
  return(res)
}

## ----lr_strata----------------------------------------------------------------
lr <- learner$new(y ~ sex + strata(trt),
                  estimate=est, predict=pred, specials = "strata")
des <- lr$design(head(pbc))
des
des$strata

## ----lr_strata_est------------------------------------------------------------
lr$estimate(pbc)
lr
lr$predict(head(pbc))

