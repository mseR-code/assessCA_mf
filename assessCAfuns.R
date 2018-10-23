# Functions for assessCA TMB model

# Plot multi-panel time series plots
plotMulti <- function(  ts = c("Nt","Bt","Ft"),
                        report = repRE, initYear = 1965 )
{
  par(mfrow = c(length(ts),1), oma = c(3,4,1,1), mar = c(1,1,1,1) )

  argList <- list(  report = report, initYear = initYear, noPar = TRUE )

  for( tsName in ts )
  {
    plotCall <- paste("plot",tsName,sep = "")
    do.call(  what = eval(plotCall), 
              args = argList, 
              envir = .GlobalEnv )
  }
}

# Plot biomass
plotSBt <- function( report = repFE, initYear = 1965, noPar = FALSE )
{
  # Pull stuff from report
  SBt      <- report$SB_t
  
  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(SBt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    axis(side = 2, las =1 )
    axis(side = 1)
    box()
    # Plot recruitment
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    lines( x = years[1:nT], y = SBt[1:nT], lwd = 2, col = "red" )
    mtext(side = 2, text = "Biomass (kt)", line = 3)
}

# Plot fishing mortality
plotNt <- function( report = repFE, initYear = 1965, noPar = FALSE )
{
  # Pull stuff from report
  Nat       <- report$N_at
  Nt        <- colSums(Nat, na.rm = T)
  
  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(Nt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    axis(side = 2, las =1 )
    axis(side = 1)
    box()
    # Plot recruitment
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    lines( x = years, y = Nt, lwd = 2, col = "grey40" )
    mtext(side = 2, text = "Numbers (1e6)", line = 3)
}


# Plot fishing mortality
plotFtg <- function( report = repFE, initYear = 1965, noPar = FALSE )
{
  # Pull stuff from report
  Ftg     <- report$F_tg

  nG      <- report$nG
  nT      <- report$nT

  cols    <- brewer.pal( nG, "Dark2" )

  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot F series
  plot( x = range(years), y = c(0,max(Ftg, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    axis(side = 2, las =1 )
    axis(side = 1)
    box()
    # Plot recruitment
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    for(g in 1:nG)
      lines( x = years[1:nT], y = Ftg[,g], lwd = 2, col = cols[g] )
    
    mtext(side = 2, text = "Fishing Mortality (/yr)", line = 3)
}

plotSag <- function( report = repFE )
{
  # pull selectivity from report
  Sag     <- report$sel_ag

  # model dimensions
  nA      <- report$nA
  nG      <- report$nG

  # 
}

# Plot stock-recruitment curve
plotSR <- function( report = repFE )
{
  # Pull stuff from report
  Rt      <- report$R_t
  Bt      <- report$SB_t
  omegaRt <- report$omegaR_t
  sigmaR  <- report$sigmaR
  R0      <- report$R0
  B0      <- report$B0
  rec.a   <- report$rec_a
  rec.b   <- report$rec_b
  h       <- round(report$rSteepness,2)
  phi     <- report$phi
  
  # Get number of time steps
  nT      <- report$nT

  SB <- seq(0,B0,length = 1000 )
  R  <- rec.a * SB / (1 + rec.b*SB)

  B20 <- 0.2*B0
  R20 <- rec.a * B20 / (1 + rec.b*B20)

  plot( x = range(SB), y = range(R,Rt,na.rm = T ),
        type = "n", las = 1, xlab = "",
        ylab = "")
    mtext( side = 1, text = "Spawning Biomass (kt)", line = 2.5)
    mtext( side = 2, text = "Recruits (1e6)", line = 2.5)
    lines( x = SB, y = R, lwd = 3 )
    points( x = Bt[1:nT], y = Rt[2:(nT+1)], pch = 16, col = "grey60" )
    # Plot B0,R0
    segments( x0 = B0, x1 = B0, y0 = 0, y1 = R0,
              lty = 2, lwd = 2 )
    segments( x0 = 0, x1 = B0, y0 = R0, y1 = R0,
              lty = 2, lwd = 2 )
    # Plot B20,R20
    segments( x0 = B20, x1 = B20, y0 = 0, y1 = R20,
              lty = 2, lwd = 2 )
    segments( x0 = 0, x1 = B20, y0 = R20, y1 = R20,
              lty = 2, lwd = 2 )
    # Label with steepness
    panLab( x = 0.8, y = 0.95, txt = paste("h = ", h, sep = "") )

}

# recruitments, adjusted for brood year
plotRt <- function( report = repFE, initYear = 1965, noPar = FALSE )
{
  # Pull stuff from report
  Rt      <- report$R_t
  omegaRt <- report$omegaR_t
  sigmaR  <- report$sigmaR
  R0      <- report$R0

  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)

  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(Rt, na.rm =T) ),
        type = "n", xlab = "", ylab = "Recruitments (1e6)",
        las = 1, axes = FALSE )
    axis(side = 2, las =1 )
    axis(side = 3)
    box()
    # Plot recruitment
    abline( h = R0, lty = 2, lwd = 1)
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    lines( x = years, y = Rt, lwd = 2, col = "grey40" )
    points( x = years, y = Rt, lwd = 2, col = "grey40" )
    mtext(side = 2, text = "Recruitments (1e6)", line = 3)

}

# recruitments, adjusted for brood year
plotRtResids <- function( report = repFE, initYear = 1965, noPar = FALSE )
{
  # Pull stuff from report
  Rt      <- report$R_t
  omegaRt <- report$omegaR_t
  sigmaR  <- report$sigmaR
  R0      <- report$R0

  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)

  # Adjust by brood year
  vertLines <- seq(from = initYear, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitment resids
  plot( x = range(years), y = range(omegaRt,omegaRt-sigmaR^2/2, na.rm =T),
        type = "n", xlab = "", ylab = "Recruitment log-residuals",
        las = 1 )
    # Plot recruitment
    abline( h = 0, lty = 2, lwd = 1)
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    abline( h = mean(omegaRt,na.rm =T), lty = 3, lwd = 2, col = "red")
    points( x = years[1:(nT)], y = omegaRt, col = "grey40", pch = 16 )
    mtext(side = 2, text = "Recruitment log-resids", line = 3)  
    panLegend(  x = 0.1, y = 0.95, 
                bty = "n",
                legTxt = c("Mean Resid"),
                lty = c(3),
                lwd = c(2),
                col = c("red") )
}

# plotBtDat <- function(  report = repFE, initYear = 1965, noPar = FALSE,
#                         surveyNames = c("Trap","Std.","StRS") )
# {
#   # Pull stuff from report
#   Bt      <- report$B_t
#   Ct      <- report$C_t
#   Itg     <- report$I_tg
#   qg      <- report$qhat_g
#   nG      <- report$nG

#   # Scale indices by q
#   ItgScaled <- Itg
#   for( g in 1:nG)
#     ItgScaled[,g] <- Itg[,g] / qg[g]

#   ItgScaled[ItgScaled < 0] <- NA

#   # Adjust by brood year
#   kage    <- report$kage

#   # Get number of time steps
#   nT    <- length(Bt)
#   years <- seq(from = initYear, length = nT, by = 1)

#   vertLines <- seq(from = initYear, to = max(years), by = 10)
  
#   # Colours for survey data
#   cols <- brewer.pal(nG, name = "Dark2")

#   # Set up plotting area
#   if(!noPar)
#     par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

#   # Plot BtCtIt
#   plot( x = range(years), y = c(0,max(Bt,ItgScaled, na.rm =T) ),
#         type = "n", xlab = "", ylab = "Biomass (kt)",
#         las = 1 )
#     abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
#     # Plot catch
#     rect( xleft = years - .3, xright = years + .3,
#           ybottom = 0, ytop = Ct, col = "grey60", border = NA )
#     # Plot biomass
#     lines( x = years, y = Bt, lwd = 3, col = "red" )
#     # Plot scaled indices
#     for( g in 1:nG)
#       points( x = years[-nT], y = ItgScaled[,g], pch = g + 20, 
#               cex = 1, bg = cols[g] )
#     panLegend(  x = 0.8, y = .95,
#                 legTxt = surveyNames,
#                 pch = 20 + 1:nG,
#                 pt.bg = cols, bty = "n" )
#     mtext(side = 2, text = "Biomass (kt)", line = 3)

# }


# Plot Biomass, catch and indices on the same
# axes
plotItResids <- function( report = repFE, initYear = 1965, noPar = FALSE,
                          surveyNames = c("Trap", "Std.", "StRS") )
{
  # Pull stuff from report
  Bt  <- report$B_t
  Itg <- report$I_tg
  qg  <- report$qhat_g
  nG  <- report$nG

  # Scale indices by q
  ItgScaled <- Itg
  for( g in 1:nG)
    ItgScaled[,g] <- Itg[,g] / qg[g]

  ItgScaled[ItgScaled < 0] <- NA

  cols <- brewer.pal(nG, name = "Dark2")

  # Get number of time steps
  nT    <- length(Bt)
  years <- seq(from = initYear, length = nT, by = 1)

  resids <- ItgScaled

  # Calc residuals
  for( g in 1:nG )    
    resids[,g] <- log(ItgScaled[,g] / Bt[-nT] )

  # Now set up plotting window
  if( !noPar )
    par( mfrow = c(1, 1), mar = c(2,3,1,3), oma = c(3,3,1,1) )

  plot( x = range(years), y = range(resids,na.rm=T), xlab = "",
          ylab = "log-Residual", type = "n", las = 1 )
    # show 0 line
    abline( h = 0, lty = 3, lwd = .8 )
    # Now show the resids
    for(g in 1:nG)
      points( x = years[-nT], y = resids[,g], pch = g + 20,
              col = cols[g] )
    panLegend(  x = 0.2, y = .95,
                legTxt = surveyNames,
                pch = 20 + 1:nG,
                col = cols, bty = "n" )
    mtext(  side = 2, outer = F,
            text = paste("Survey log-residuals", sep = ""),
            line = 3 )
}

