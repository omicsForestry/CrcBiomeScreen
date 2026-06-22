#' Internal caret model specification for XGBoost
#' This internal helper defines a caret-compatible model specification for
#' XGBoost. Recent versions of xgboost return booster objects containing
#' `XGBAltrepPointerClass` objects. These objects may fail when
#' `caret::train()` attempts to append metadata, such as `xNames`, directly
#' to the fitted model object.
#' This helper preserves the caret training, resampling, and tuning workflow,
#' but wraps the fitted xgboost booster inside a standard R list. It also
#' defines compatible `predict()` and `prob()` methods so that caret can
#' obtain class predictions and class probabilities without directly modifying
#' the xgboost booster object.
#' The function is used internally by the XGBoost modelling functions and is
#' not intended to be called directly by users.
#' @return A list defining a custom caret model specification for XGBoost.
#' @keywords internal
#' @noRd

.getCaretXgbTreeCompat <- function() {
  custom_xgb <- caret::getModelInfo("xgbTree", regex = FALSE)[[1]]
  custom_xgb$label <- "Compatible XGBoost Tree"

  custom_xgb$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    if (!inherits(x, "xgb.DMatrix")) {
      x <- as.matrix(x)
    }

    if (is.factor(y)) {
      if (length(lev) == 2) {
        y <- ifelse(y == lev[1], 1, 0)

        if (!inherits(x, "xgb.DMatrix")) {
          x <- xgboost::xgb.DMatrix(x, label = y, missing = NA)
        } else {
          xgboost::setinfo(x, "label", y)
        }

        if (!is.null(wts)) {
          xgboost::setinfo(x, "weight", wts)
        }

        out <- xgboost::xgb.train(
          params = list(
            eta = param$eta,
            max_depth = param$max_depth,
            gamma = param$gamma,
            colsample_bytree = param$colsample_bytree,
            min_child_weight = param$min_child_weight,
            subsample = param$subsample,
            objective = "binary:logistic"
          ),
          data = x,
          nrounds = param$nrounds,
          verbose = 0,
          ...
        )
      } else {
        y <- as.numeric(y) - 1

        if (!inherits(x, "xgb.DMatrix")) {
          x <- xgboost::xgb.DMatrix(x, label = y, missing = NA)
        } else {
          xgboost::setinfo(x, "label", y)
        }

        if (!is.null(wts)) {
          xgboost::setinfo(x, "weight", wts)
        }

        out <- xgboost::xgb.train(
          params = list(
            eta = param$eta,
            max_depth = param$max_depth,
            gamma = param$gamma,
            colsample_bytree = param$colsample_bytree,
            min_child_weight = param$min_child_weight,
            subsample = param$subsample,
            objective = "multi:softprob",
            num_class = length(lev)
          ),
          data = x,
          nrounds = param$nrounds,
          verbose = 0,
          ...
        )
      }
    } else {
      if (!inherits(x, "xgb.DMatrix")) {
        x <- xgboost::xgb.DMatrix(x, label = y, missing = NA)
      } else {
        xgboost::setinfo(x, "label", y)
      }

      if (!is.null(wts)) {
        xgboost::setinfo(x, "weight", wts)
      }

      out <- xgboost::xgb.train(
        params = list(
          eta = param$eta,
          max_depth = param$max_depth,
          gamma = param$gamma,
          colsample_bytree = param$colsample_bytree,
          min_child_weight = param$min_child_weight,
          subsample = param$subsample,
          objective = "reg:squarederror"
        ),
        data = x,
        nrounds = param$nrounds,
        verbose = 0,
        ...
      )
    }

    list(
      model = out,
      obsLevels = lev,
      lev = lev,
      problemType = if (is.factor(y)) "Classification" else "Regression"
    )
  }

  custom_xgb$predict <- function(modelFit, newdata, submodels = NULL) {
    if (!inherits(newdata, "xgb.DMatrix")) {
      newdata <- as.matrix(newdata)
      newdata <- xgboost::xgb.DMatrix(newdata, missing = NA)
    }

    pred <- predict(modelFit$model, newdata)

    lev <- modelFit$obsLevels
    if (is.null(lev)) {
      lev <- modelFit$lev
    }

    if (!is.null(lev) && length(lev) == 2) {
      out <- ifelse(pred >= 0.5, lev[1], lev[2])
      return(factor(out, levels = lev))
    }

    pred
  }

  custom_xgb$prob <- function(modelFit, newdata, submodels = NULL) {
    if (!inherits(newdata, "xgb.DMatrix")) {
      newdata <- as.matrix(newdata)
      newdata <- xgboost::xgb.DMatrix(newdata, missing = NA)
    }

    pred <- predict(modelFit$model, newdata)

    lev <- modelFit$obsLevels
    if (is.null(lev)) {
      lev <- modelFit$lev
    }

    if (!is.null(lev) && length(lev) == 2) {
      out <- data.frame(
        class1 = pred,
        class2 = 1 - pred,
        check.names = FALSE
      )
      colnames(out) <- lev
      return(out)
    }

    stop("Probability prediction for multiclass XGBoost is not yet implemented.")
  }

  custom_xgb$loop <- NULL

  custom_xgb
}
