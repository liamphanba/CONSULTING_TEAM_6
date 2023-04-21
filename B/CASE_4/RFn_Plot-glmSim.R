plot.glmSim <- function(obj, which=c(1L:3L), rob=FALSE, SEED=NULL, Nsim=19,
                        smoother=function(x,y) lowess(x, y, f=2/3, iter=0),
                        singular.ok=TRUE)
{
    ## rkst, Version Oct/12/2021
    if(!inherits(obj, "glm")) stop("use only with \"glm\" objects")
    show <- rep(FALSE, 3)
    show[which] <- TRUE

    X <- model.matrix(obj)
    x.n <- nrow(X)
    etaH <- predict(obj, type="link")
    ## w <- weights(obj)
    ## if(is.null(w)) w <- rep(1, nrow(X))
    pWeights <- weights(obj, type="prior")
    binaryResp <- (family(obj)$family == "binomial") & any(pWeights==1)

    phi <- summary(obj)$dispersion
    OffSet <- obj$offset
    Family <- family(obj)
    Control <- obj$control
    muH <- predict(obj, type="response")

    y.sim <- switch(obj$family$family,
                    "binomial" = {function(){
                        rbinom(x.n,size=pWeights, prob=muH)/pWeights}},
                    "gaussian" = function() rnorm(x.n,muH,sqrt(phi)),
                    "Gamma" = {function()
                        rgamma(x.n, shape=1/phi, scale= muH*phi)},
                    "inverse.gaussian" = {
                        require(statmod)
                        warning("under construction")
                        function() rinvgauss(x.n, mean=muH,
                                                dispersion=phi)},
                    "poisson" = function() rpois(x.n,lambda=muH),
                    stop("not (yet) implemented"))

    ## Tukey-Anscombe Plot
    if(show[1L]){
        res <- residuals.glm(obj, type="pearson")
        ylim <- extendrange(r=range(res, na.rm = TRUE), f = 0.08)
        plot(etaH, res, xlab="Fitted Linear Predictor",
             ylab = "Pearson Residuals", main="", ylim=ylim, type="n")
        abline(h=0, lty=3) ## , col = "gray"
        lines(smoother(etaH, res), lwd=1.5, col="red")

        if(!is.null(SEED)) set.seed(SEED)
        for(i in 1:Nsim){
            FIT <- glm.fit(x=X, y=y.sim(), weights=pWeights, start=coef(obj),
                           offset=OffSet, family=Family, control=Control,
                           intercept=FALSE, singular.ok=singular.ok)
            ## class(FIT) <- c(FIT$class, c("glm", "lm"))
            ## lines(lowess(predict(FIT), resid(FIT), f=2/3, iter=3),col="grey")
            lines(smoother(predict.glm(FIT),
                           residuals.glm(FIT, type="pearson")), col="grey")
        }
         lines(smoother(etaH, res), lwd=1.5, col="red")
   }
    if(any(show[2L:3L])){
         f.stdres <- function(glmFit, pW=pWeights) {
            r <- residuals.glm(glmFit, type="pearson")
            hPhi <- summary.glm(glmFit)$dispersion
            hHii <- lm.influence(glmFit, do.coef=FALSE)$hat
            ## Be careful, there are problems when hii==1
            sr <- r*sqrt(pW/(hPhi*(1 - hHii)))
            sr[hHii >= 1] <- NaN
            sr
        }
    }

    ## normal plot
    if(show[2L]){
        SIM <- matrix(NaN, ncol=Nsim, nrow=x.n)
        if(Nsim>0){
            if(!is.null(SEED)) set.seed(SEED)
            for(i in 1: Nsim){
                FIT <- glm.fit(x=X, y=y.sim(), weights=pWeights,
                               start=coef(obj), offset=OffSet, family=Family,
                               control=Control, intercept=FALSE,
                               singular.ok=singular.ok)
                SIM[,i] <- sort(f.stdres(FIT))
            }
        }
        RQQN <- qqnorm(f.stdres(obj), plot.it=FALSE)
        ylim <- range(c(RQQN$y, SIM), finite=TRUE)
        plot(range(RQQN$x, finite=TRUE), ylim, type="n",
             xlab = "Theoretical Quantiles", ylab = "Std. Pearson Resid.")
        if(Nsim>0)
            points(rep(sort(RQQN$x),Nsim), as.vector(SIM), col="gray")
        points(RQQN$x, RQQN$y, lwd=2)
    }
    ## Scale-Location Plot
    if(show[3L]){
        sqrtabsR <- sqrt(abs(f.stdres(obj)))
        ok <- is.finite(sqrtabsR)
        ylim <- c(0, max(sqrtabsR[ok]))
        yl <- as.expression(substitute(sqrt(abs(YL)),
                                 list(YL=as.name("Std. Pearson Resid."))))
        plot(etaH, sqrtabsR, xlab="Fitted Linear Predictor", ylab=yl, main="",
             ylim=ylim, type="n")
        lines(smoother(etaH[ok], sqrtabsR[ok]), lwd=1.5, col="red")
        if(!is.null(SEED)) set.seed(SEED)
        for(i in 1:Nsim){
            FIT <- glm.fit(x=X, y=y.sim(), weights=pWeights,
                           start=coef(obj), offset=OffSet, family=Family,
                           control=Control, intercept=FALSE,
                           singular.ok=singular.ok)
            class(FIT) <- c(FIT$class, c("glm", "lm"))

            h.sasRes  <- sqrt(abs(f.stdres(FIT)))
            ok <- is.finite(h.sasRes)
            lines(smoother(predict(FIT)[ok], h.sasRes[ok]), col="grey")
        }
         lines(smoother(etaH[ok], sqrtabsR[ok]), lwd=1.5, col="red")
   }
    invisible()
}
