
genDynamicOccData <- function(nsite = 100,
                              nreps = 5,
                              nyear = 15,
                              mu.p,
                              psi1,
                              sigma.p,
                              mu.phi,
                              sigma.phi,
                              mu.gamma,
                              sigma.gamma) {

    ##  Generation and analysis of simulated data for multi season
    ##  occupancy model (adapted from Kery and Schaud 2012)
    ## Set up some required arrays
    site <- 1:nsite					## Sites
    year <- 1:nyear					## Years
    psi <- rep(NA, nyear)				## Occupancy probability
    muZ <- z <- array(dim = c(nsite, nyear))	## Expected and realized occurrence
    y <- array(NA, dim = c(nsite, nreps, nyear))	## Detection histories
    ## Determine initial occupancy and demographic parameters
    psi[1] <- psi1				## Initial occupancy probability
    p <-  rnorm(nyear, mu.p, sigma.p)
    phi <- rnorm(nyear -1, mu.phi, sigma.phi)
    gamma <- rnorm(nyear -1, mu.gamma, sigma.gamma)
    ## Generate latent states of occurrence
    ## First year
    z[,1] <- rbinom(nsite, 1, psi[1])		## Initial occupancy state
    ## Later years
    for(site in 1:nsite){				## Loop over sites
        for(year in 2:nyear){				## Loop over years
            muZ[site, year] <- z[site, year-1]*expit(phi[year-1]) +
                (1-z[site, year-1])*expit(gamma[year-1]) ## Prob for occ.
            z[site,year] <- rbinom(1, 1, muZ[site, year])
        }
    }
    ## Generate detection/nondetection data
    for(site in 1:nsite){
        for(year in 1:nyear){
            prob <- z[site, year]*expit(p[year])
            for(rep in 1:nreps){
                y[site, rep ,year] <- rbinom(1, 1, prob)
            }
        }
    }

    ## Compute annual population occupancy
    for (year in 2:nyear){
        psi[year] <- psi[year-1]*expit(phi[year-1]) +
            (1-psi[year-1])*expit(gamma[year-1])
    }

    return(list(nsite = nsite,
                nreps = nreps,
                nyear = nyear,
                psi = psi,
                z = z,
                phi = phi,
                gamma = gamma,
                p = p,
                y = y,
                mu.p = mu.p,
                sigma.p = sigma.p,
                mu.phi = mu.phi,
                sigma.phi = sigma.phi,
                mu.gamma = mu.gamma,
                sigma.gamma = sigma.gamma))
}
