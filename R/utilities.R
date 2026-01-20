#' @rdname utilities
#' @title Utilities
#' @name step_s
#' @description Given a data-frame with an "s" step function, expand the data-frame
#' to include the steps.
#' @param data data-frame
#' @param x name of the x-variable (required)
#' @param y name of the y-variable (required)
#' @param ymin name of the ymin variable (optional)
#' @param ymax name of the ymax variable (optional)
#' @param group name of a grouping variable (optional)
#' @param add_origin logical for whether to add an origin to the start of each group
#' @param x_origin double for the value of x at the origin
#' @param y_origin double for the value of y, ymin and ymax at the origin
#' @return expanded data-frame with the same names
#' @importFrom utils head
#' @export
#' @examples
#'   step_s(data.frame(g=c(1,1), a=1:2, b=4:5), a, b)
#'   step_s(data.frame(g=c(2,2), a=3:4, b=6:7), a, b)
#'   step_s(data.frame(g=c(1,1,2,2), a=1:4, b=4:7), a, b, group=g)
step_s = function(data,x,y,ymin,ymax,group,add_origin=TRUE,x_origin=0,y_origin=1) {
    data = as.data.frame(data)
    .call = match.call()
    if (!missing(group)) {
        nm = deparse(substitute(group))
        newdata = by(data, data[[nm]], function(datai) {
            .call$data = datai
            .call$group=NULL
            eval(.call)
        })
        newdata = do.call(rbind, newdata)
        rownames(newdata) = NULL
        newdata
    } else {
        xn = deparse(substitute(x))
        yn = deparse(substitute(y))
        index = 1:nrow(data)
        if (add_origin) {
            xindex=c(1,rep(index,each=2))
            yindex=head(c(1,1,rep(index,each=2)),-1)
        } else {
            xindex=rep(index,each=2)
            yindex=head(rep(index,each=2),-1)
        }
        newdata = data[xindex,,drop=FALSE]
        newdata[[yn]] = data[[yn]][yindex]
        if (add_origin) {
            newdata[[xn]][1] = x_origin
            newdata[[yn]][1:2] = y_origin
        }
        if (!missing(ymin)) {
            nm = deparse(substitute(ymin))
            newdata[[nm]] = data[[nm]][yindex]
            if (add_origin)
                newdata[[nm]][1:2] = y_origin
        }        
        if (!missing(ymax)) {
            nm = deparse(substitute(ymax))
            newdata[[nm]] = data[[nm]][yindex]
            if (add_origin)
                newdata[[nm]][1:2] = y_origin
        }        
        rownames(newdata) = NULL
        newdata
    }
}

#' @rdname utilities
#' @name as.data.frame.survfit
#' @description Given a survfit object, return a data-frame
#' @param x a survfit object
#' @param row.names not used (in generic signature)
#' @param optional not used (in generic signature)
#' @param type which type of data-frame to return: either expanded (e.g. for plotting) or plain (not expanded)
#' @param ... not used (in generic signature)
#' @export
#' @examples
#'   library(survival)
#'   library(tinyplot)
#'   sfit1 = survfit(Surv(time,status)~rx, data=survival::colon,
#'                   subset=etype==1)
#'   with(as.data.frame(sfit1),
#'        tinyplot::plt(surv~time|strata,ymin=lower,ymax=upper,type="ribbon"))
as.data.frame.survfit = function(x, row.names, optional, type=c("expanded","plain"), ...) {
    type = match.arg(type)
    strata <- time <- surv <- lower <- upper <- NULL
    if ("strata" %in% names(x)) {
        df = with(x, data.frame(strata=rep(names(strata), strata),
                                time, surv, lower, upper))
        if (type=="plain") df else step_s(df, time, surv, lower, upper, group=strata)
    } else {
        df = with(x, data.frame(time, surv, lower, upper))
        if (type=="plain") df else step_s(df, time, surv, lower, upper)
    }
}

#' @rdname utilities
#' @name as.data.frame.summary.survfit
#' @description Given a summary.survfit object, return a data-frame
#' @param x a summary.survfit object
#' @param row.names not used (in generic signature)
#' @param optional not used (in generic signature)
#' @param type which type of data-frame to return: either expanded (e.g. for plotting) or plain (not expanded)
#' @param ... not used (in generic signature)
#' @export
#' @importFrom stats qnorm
#' @examples
#'   library(survival)
#'   library(tinyplot)
#'   sfit1 = survfit(Surv(time,status)~rx, data=survival::colon,
#'                   subset=etype==1)
#'   with(as.data.frame(sfit1, type="expanded"),
#'        tinyplot::plt(surv~time|strata,ymin=lower,ymax=upper,type="ribbon"))
as.data.frame.summary.survfit = function(x, row.names, optional,
                                         type=c("expanded","plain"), ...) {
    strata <- time <- surv <- lower <- upper <- NULL
    type = match.arg(type)
    if ("strata" %in% names(x)) {
        df = with(x, data.frame(strata=rep(names(strata), strata),
                                time, surv, lower, upper))
        if (type=="plain") df else step_s(df, time, surv, lower, upper, group=strata)
    } else {
        df = with(x, data.frame(time, surv, lower, upper))
        if (type=="plain") df else step_s(df, time, surv, lower, upper)
    }
}

#' @rdname utilities
#' @name predict_coxph_tt
#' @title Calculate the linear predictor for the tt argument for a
#'     coxph fit
#' @description Calculates lp=tt(x). Typically, the tt function should
#'     include an intercept term (see the examples below). Note that
#'     spline terms assume that the x argument is multiplicative;
#'     moreover, the additional arguments are not passed. For other
#'     types of tt terms, the x is passed directly to the tt function
#'     together with other arguments.
#' @param object a coxph fit with one (and currently only one) tt
#'     argument
#' @param times a numeric vector of times to evaluate the linear
#'     predictor
#' @param type a character for the type of prediction (currently only
#'     the linear predictor for the tt argument)
#' @param se.fit a logical for whether to return the standard errors
#' @param x a numeric scalar to pass to the tt function (assumed to be
#'     multiplicative for any spline terms)
#' @param ... other parameters to pass to the tt function
#' @return a vector of fitted values (when se.fit=FALSE) or a data-frame
#'     with fitted and se.fit columns (when se.fit=TRUE)
#' @examples
#'  library(splines)
#'  library(tinyplot)
#'  fit1 = coxph(Surv(time,status)~tt(treat),data=breast_rfs,
#'               tt=function(x,t,...) x*cbind(1,t))
#'  fit2 = coxph(Surv(time,status)~tt(treat),data=breast_rfs,
#'               tt=function(x,t,...) x*ns(t,df=4,intercept=TRUE))
#'  times = seq(0,2500,len=301L)
#'  df1 = transform(predict_coxph_tt(fit1,times,se.fit=TRUE),
#'                  lower=exp(fitted-1.96*se.fit),
#'                  upper=exp(fitted+1.96*se.fit),
#'                  fitted=exp(fitted),
#'                  model="linear",times=times)
#'  df2 = transform(predict_coxph_tt(fit2,times,se.fit=TRUE),
#'                  lower=exp(fitted-1.96*se.fit),
#'                  upper=exp(fitted+1.96*se.fit),
#'                  fitted=exp(fitted),
#'                  model="ns",times=times)
#'  with(rbind(df1,df2),
#'       plt(fitted ~ times | model, ymin=lower, ymax=upper, type="ribbon",
#'           xlab="Time since diagnosis (days)",
#'           ylab="Hazard ratio comparing treated with untreated"))
#'  with(subset(breast_rfs,status==1), rug(time))
#' @importFrom stats vcov coef
#' @export
predict_coxph_tt = function(object, times, type="lp", se.fit=FALSE, x=1, ...) {
    type = match.arg(type)
    tt_index = attr(object$terms,"specials")$tt
    stopifnot(inherits(object,"coxph"),
              is.integer(tt_index),
              length(tt_index)==1)
    tt_predvars = attr(object$terms,"predvars")[[tt_index+1]]
    if (tt_predvars[[1]] == as.name("tt")) {
        X = eval(object$call$tt)(x,times,...)
    } else {
        assign(as.character(tt_predvars[[2]]), times)
        X = x*eval(tt_predvars)
    }
    if (is.null(attr(X,"dim")))
        X = matrix(X,ncol=1)
    n = length(coef(object))
    index = (n-ncol(X)+1):n
    tt_coef = coef(object)[index]
    fitted = drop(X %*% tt_coef)
    if (!se.fit)
        return(fitted)
    ## check = sqrt(diag(X %*% vcov(object)[index,index] %*% t(X)))
    se.fit = sqrt(rowSums(X %*% vcov(object)[index,index] * X))
    data.frame(fitted,se.fit)
}

#' @rdname utilities
#' @name predict_coxph_tv
#' @title Predicted survival from a coxph object for subjects with time-varying exposures
#' @param object a coxph object
#' @param data a data-frame or tmerge object with columns tstart,
#'     tstop, status, the exposure variables and an id variable
#' @param id a character for the subject id
#' @returns an update of data with survival probabilities
#' @examples
#' library(survival)
#' liver = transform(collett::liverbase, lbr=NULL)
#' liver = tmerge(liver, liver, id=patient, status=event(time,status))
#' liver = tmerge(liver,
#'                rbind(with(collett::liverbase, data.frame(patient,tstart=0,lbr)),
#'                      with(collett::lbrdata0, data.frame(patient,tstart=time,lbr))),
#'            id=patient, lbr = tdc(tstart,lbr))
#' fit3 = coxph(Surv(tstart,tstop,status)~lbr+treat,liver)
#' predict_coxph_tv(fit3,data=subset(liver,patient %in% c(1,7)),
#'                  id="patient")
#' @export 
predict_coxph_tv = function(object,data,id) {
    stopifnot(inherits(object,"coxph"),
              inherits(data,"data.frame"),
              all(c("tstart","tstop","status",id) %in% names(data)))
    data = data[order(cbind(data[[id]],data$tstart)),]
    data_list = by(data, data[[id]], function(datai) {
        ## interval-specific survival
        interval_S = by(datai, 1:nrow(datai),
                        function(newdata) {
                            times = c(newdata$tstart, newdata$tstop)
                            surv = summary(survfit(object, newdata=newdata),
                                           times=times)$surv
                            surv[2]/surv[1]
                        })
        datai$surv = cumprod(interval_S)
        datai
    })
    do.call(rbind, data_list)
}


#' @rdname utilities
#' @name plot_coxph_functional
#' @title Plot the functional form for a coxph model using smoothed martingale residuals
#' @param formula a formula with a Surv on the lhs and a single variable on the rhs
#' @param data a dataset for evaluation of the coxph model
#' @param x a numeric vector for the smoother (defaults to the 401 values between the range)
#' @param pch an integer for the pch argument in the plot for the residuals
#' @param ylab a character for the ylab argument in the plot
#' @param smoother a character for the name of the smoother
#' @param smoother.formula a formula for the smoother in terms of resi and xi
#' @param smoother.args a list of arguments to pass to the smoother function
#' @param points.args a list of arguments to pass to the points function
#' @param ... other arguments to pass to the plot function
#' @returns invisible plot return
#' @import survival
#' @import tinyplot
#' @importFrom stats residuals predict qnorm
#' @importFrom utils modifyList
#' @importFrom survival coxph
#' @importFrom graphics points
#' @export
#' @examples
#' library(survival)
#' par(mfrow=c(2,2))
#' plot_coxph_functional(Surv(time,status)~hb, data=collett::myeloma,
#'                       xlab="Value of Hb")
#' plot_coxph_functional(Surv(time,status)~bun, data=collett::myeloma,
#'                       xlab="Value of Bun")
#' plot_coxph_functional(Surv(time,status)~log(bun), data=collett::myeloma,
#'                       xlab="Value of log Bun")
plot_coxph_functional = function(formula,
                            data,
                            x=NULL,
                            pch=19,
                            ylab="Martingale residual for null model",
                            smoother=c("loess","lm"),
                            smoother.formula = resi ~ xi,
                            smoother.args=list(),
                            points.args=list(),
                            ...) {
    smoother = match.arg(smoother)
    smoother.formula = substitute(smoother.formula)
    stopifnot(length(formula)==3,
              is.null(x) || is.numeric(x))
    xi = eval(formula[[3]], data)
    formula[[3]] = 1
    resi = residuals(coxph(formula, data), type="martingale")
    fit1 = do.call(smoother, modifyList(list(formula=smoother.formula), smoother.args))
    if (is.null(x)) {
        x = seq(min(xi),max(xi),length=401L)
    }
    qq=qnorm(0.975)
    newdata = data.frame(xi=x)
    plt1 = with(switch(smoother,
                        loess=predict(fit1, newdata=newdata, se=TRUE),
                        lm=predict(fit1,  newdata=newdata, se.fit=TRUE)),
                plt(fit~x, ymin=fit-qq*se.fit, ymax=fit+qq*se.fit,
                    type="ribbon", ylab=ylab,
                    ylim=range(resi), ...))
    do.call(points, modifyList(list(xi, resi, pch=pch), points.args))
    invisible(plt1)
}
