#panel.lm
panel.lm <- 
function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.lm = "blue",  ...) 
{   ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin,xmin),max(ymax,xmax))
    xlim <- ylim
    points(x, y, pch = pch, col = col, bg = bg, cex = cex,ylim = ylim, xlim= xlim)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        abline(lm(y[ok]~ x[ok]), 
            col = col.lm, ...)
}

