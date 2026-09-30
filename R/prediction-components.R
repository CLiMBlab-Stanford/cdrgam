#' Predict response and link components
#'
#' `predict_components()` returns a stable representation of response- and
#' link-scale predictions for artifact writers and other downstream tools.
#' Ordinary models return vectors. Distributional models additionally return
#' matrices with one column per response-distribution parameter.
#'
#' @param object A fitted `cdrgam` model.
#' @param newdata Prediction data accepted by [predict.cdrgam()].
#' @return A named list containing primary response- and link-scale estimates
#'   and standard errors. Distributional fits also include parameter names and
#'   the complete parameter matrices on both scales.
#' @export
predict_components <- function(object, newdata) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    link_prediction <- stats::predict(
        object, newdata=newdata, type='link', se.fit=TRUE
    )
    if (is.list(link_prediction) &&
            all(c('fit', 'se.fit') %in% names(link_prediction))) {
        link_value <- link_prediction$fit
        link_se <- link_prediction$se.fit
    } else {
        link_value <- link_prediction
        link_se <- rep.int(NA_real_, length(link_value))
    }
    if (!isTRUE(object$cdrgam$distributional)) {
        response_value <- object$family$linkinv(link_value)
        response_se <- link_se * abs(object$family$mu.eta(link_value))
        return(list(
            distributional=FALSE,
            prediction=response_value,
            prediction_se=response_se,
            link_prediction=link_value,
            link_prediction_se=link_se
        ))
    }
    response_prediction <- stats::predict(
        object, newdata=newdata, type='response', se.fit=TRUE
    )
    response_value <- response_prediction$fit
    response_se <- response_prediction$se.fit
    parameter_names <- object$cdrgam$parameter_names
    colnames(link_value) <- parameter_names
    colnames(link_se) <- parameter_names
    colnames(response_value) <- parameter_names
    colnames(response_se) <- parameter_names
    list(
        distributional=TRUE,
        parameter_names=parameter_names,
        prediction=response_value[, 'location'],
        prediction_se=response_se[, 'location'],
        link_prediction=link_value[, 'location'],
        link_prediction_se=link_se[, 'location'],
        response_parameters=response_value,
        response_parameter_se=response_se,
        link_parameters=link_value,
        link_parameter_se=link_se
    )
}
