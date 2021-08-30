        # A script to simulate a stage-structured population of foothills palo verde
        # at Tumamoc Hill, with intraspecific density-dependence (alpha),
        # and then to perturb recruitment rates to represent varying possible effects of buffel grass invasion. 
        
        # =========== Estimating parameters based on empirical data available =============
        
        # s11: survival of seedlings based on data from Shreve (1917):
        total <- c(303,731,203,251,87,26,16,45) #total seedlings in each year
        survived <- c(189,81,100,53,19,16,16,19) #total of those seedlings that survived to next year
        weighted_arith_mean_age_agnostic <- round(sum(survived)/sum(total),3)
        
        # f13: Annual emergence of seedlings per adult, based on Bowers & Turner (2002) and Rodriguez-Buritica et al. (2013)
        seedlings_per_year <- c(155,163,36,170,329,155) #emergence data
        adults_area_A <- 286*(557/10000)     
        seedlings_per_adult <- round(mean(seedlings_per_year)/adults_area_A,1)
        
        # s33: Annual survival of adults based on Bowers & Turner (2002)
        
        # Enter data from B&T (2002) Fig 4 demography of age class (max age) and number of trees:
        Age <- seq(30,155,5)
        Count <- c(32, 37, 49, 50, 45, 52, 55, 50, 69, 50, 42, 34, 35, 30, 17, 13, 24, 16, 7, 13, 1, 6, 4, 1, 1, 3)
        bandt <- data.frame(Age,Count)
        
        # Smooth out stochastic variation among age classes by fitting a linear model:
        lm.bandt <- lm(bandt$Count~bandt$Age)
        smry.bt <- summary(lm.bandt)
        # and predicting counts of population at that age using the model:
        bandt$pre <- predict(lm.bandt)
        
        # Set up a for loop to calculate difference in counts a.........................nd ages
        bandt.age.difs <- c()
        bandt.pre.rats <- c()
        for(i in 1:(nrow(bandt)-2)) {
          bandt.age.difs <- c(bandt.age.difs, bandt$Age[i+1] - bandt$Age[i])
          bandt.pre.rats <- c(bandt.pre.rats, bandt$pre[i+1] / bandt$pre[i])
        } 
        # convert into annual proportional change in population size
        library(pracma)
        s.rates.bt <- nthroot(bandt.pre.rats,bandt.age.difs) 
        # Take unweighted mean
        bt.rate.uw <- round(mean(s.rates.bt),3)
        
        # g21: Mean seedlings per cohort that may have established from Shreve (1917)
        seedlings_per_year_shreve <- c(303,542,122,151,34,7) # first six years of data
        seedlings_surviving_at_end <- c(6,2,0,2,2,2) # seedlings that survived to end of study, at least 3 years
        establishment_per_cohort <- round(mean(seedlings_surviving_at_end/seedlings_per_year_shreve),4)
        
        # ================ Constructing a model based on estimates =======================
        
        # Step one: define parameters (and variability) in the absence of density-dependence (or at least with it reduced)
        s11 <- weighted_arith_mean_age_agnostic  # Annual survival of seedlings over their first 7 years, based on Shreve (1917)
        surv_age_agnostic <- survived/total # Bootstrapped 95% CI of g21
        B = 1000
        n = length(surv_age_agnostic)
        boot.samples.s11 = matrix(sample(surv_age_agnostic, size = B * n, replace = TRUE),B, n)
        boot.ci.s11 <- c(-1, 1) * sd(apply(boot.samples.s11, 1, mean)) + s11# Bootstrapped 95% CI of s11
        f13 <- 10.5 # Annual emergence of seedlings per adult, based on Bowers & Turner (2002)
        f_per_adult_per_year <- seedlings_per_year/adults_area_A 
        n = length(f_per_adult_per_year)
        boot.samples.f13 = matrix(sample(f_per_adult_per_year, size = B * n, replace = TRUE),B, n)
        boot.ci.f13 <- c(-1, 1) * sd(apply(boot.samples.f13, 1, mean)) + f13 # Bootstrapped 95% CI of f13
        s22 <- 0.99 # Annual survival of saplings based on being close to 100%
        # Just perturb s22 by 1% in either direction for lack of more precise data (not much room for perturbation being so close to 1)
        s33 <- 0.972 # Annual survival of adults based on Bowers & Turner (2002)
        n = length(s.rates.bt)
        boot.samples.s33 = matrix(sample(s.rates.bt, size = B * n, replace = TRUE),B, n)
        boot.ci.s33 <- c(-1, 1) * sd(apply(boot.samples.s33, 1, mean)) + s33# Calc 95% CI of annual survival rates estimated from Bowers & Turner (2002)
        g21 <- 0.0635 # Mean seedlings per cohort that may have established from Shreve (1917)
        est_per_co <- seedlings_surviving_at_end/seedlings_per_year_shreve # Bootstrapped 95% CI of g21
        n = length(est_per_co)
        boot.samples.g21 = matrix(sample(est_per_co, size = B * n, replace = TRUE),B, n)
        boot.ci.g21 <- c(-1, 1) * sd(apply(boot.samples.g21, 1, mean)) + g21
        g32 <- 0.0556 # 1/18th of saplings transition to adults based on lack of data
        # Perturb this by ~30% in either direction, from 1/25 to 1/11, because this value is poorly constrained
        
        # Step two: define density-independent (DI) matrix and conduct sensitivity and elasticity analyses
        DI <- matrix( 
          c(s11,0,f13,g21,s22,0,0,g32,s33), # the data elements 
          nrow=3,              # number of rows 
          ncol=3,              # number of columns 
          byrow = TRUE)        # fill matrix by rows
        
        library(popbio)
        eigen.analysis(DI) # This provides the elasticities and sensitivities in Table S2.2
        
        
        # Step three: simulate growth until it stabilizes
        # In Bowers & Turner (2002), their 1999 census had:
        # ~63 trees between 7 and 25, and
        # 751 trees >25 y/o on 6 ha, and
        # ~9 seedlings emerging per adult
        # So to express these in per-hectare terms, 
        # let's call initial densities:
        # 125 adults, 11 saplings, and 1125 seedlings
        
        N0 <- c(1125,11,125)
        N1 <- DI%*%N0
        N2 <- DI%*%N1
        
        # Density-depndent (DD) weights for age classes are based on horizontal area of crown
        c1 <- 3.14*(.01^2)    # seedlings are ~2cm diameter
        c2 <- 3.14*(.1^2)     # saplings up to 25 y/o when they start producing seeds could be 5-50 cm diameter
        c3 <- 3.14*(1.5^2)    # adults producing seeds would be typically 0.5-6 m in diameter
        cs <- c(c1,c2,c3)
        
        ### Conduct a few simulations of population densities with various values of alpha to see where the
        # population denisites stabilize (Fig 2.1).
        # For each value of alpha also perturb all vital rates up and down by 1% to get a sense of 
        # the variability as a result
        
        
        alpha.tmp <- .03   # Density dependence factor
        pops <- data.frame(N0)       # Part A of the figure: low alpha, perturbed F13
        pops.hi <- data.frame(N0)
        pops.lo <- data.frame(N0)
        for(i in 1:3000){     
          N.old.tmp <- pops[,i]
          N.agg.tmp <- cs*N.old.tmp
          N.old.tmp.hi <- pops.hi[,i]
          N.agg.tmp.hi <- cs*N.old.tmp.hi
          N.old.tmp.lo <- pops.lo[,i]
          N.agg.tmp.lo <- cs*N.old.tmp.lo
          DD.tmp <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp))),0,f13/(1+(alpha.tmp*sum(N.agg.tmp))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          DD.tmp.hi <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.hi))),0,(boot.ci.f13[2])/(1+(alpha.tmp*sum(N.agg.tmp.hi))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.hi))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows      
          DD.tmp.lo <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.lo))),0,(boot.ci.f13[1])/(1+(alpha.tmp*sum(N.agg.tmp.lo))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.lo))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          N.new.tmp <- DD.tmp%*%N.old.tmp
          N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
          N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
          pops <- data.frame(pops,N.new.tmp)
          pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
          pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
        }
        tpop_019_F13 <- data.frame(t(pops))
        colnames(tpop_019_F13) <- c("Seedlings","Saplings","Adults")
        tpop_019_F13$Year <- seq(0,3000,1)   
        tpop_019_F13_hi <- data.frame(t(pops.hi))
        colnames(tpop_019_F13_hi) <- c("Seedlings","Saplings","Adults")
        tpop_019_F13_hi$Year <- seq(0,3000,1)   
        tpop_019_F13_lo <- data.frame(t(pops.lo))
        colnames(tpop_019_F13_lo) <- c("Seedlings","Saplings","Adults")
        tpop_019_F13_lo$Year <- seq(0,3000,1)   
        
        pops <- data.frame(N0)       # Part B of the figure: low alpha, perturbed S11
        pops.hi <- data.frame(N0)
        pops.lo <- data.frame(N0)
        for(i in 1:3000){     
          N.old.tmp <- pops[,i]
          N.agg.tmp <- cs*N.old.tmp
          N.old.tmp.hi <- pops.hi[,i]
          N.agg.tmp.hi <- cs*N.old.tmp.hi
          N.old.tmp.lo <- pops.lo[,i]
          N.agg.tmp.lo <- cs*N.old.tmp.lo
          DD.tmp <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp))),0,f13/(1+(alpha.tmp*sum(N.agg.tmp))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          DD.tmp.hi <- matrix(
            c(boot.ci.s11[2]/(1+(alpha.tmp*sum(N.agg.tmp.hi))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.hi))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.hi))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows      
          DD.tmp.lo <- matrix(
            c(boot.ci.s11[1]/(1+(alpha.tmp*sum(N.agg.tmp.lo))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.lo))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.lo))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          N.new.tmp <- DD.tmp%*%N.old.tmp
          N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
          N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
          pops <- data.frame(pops,N.new.tmp)
          pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
          pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
        }
        tpop_019_S11 <- data.frame(t(pops))
        colnames(tpop_019_S11) <- c("Seedlings","Saplings","Adults")
        tpop_019_S11$Year <- seq(0,3000,1)   
        tpop_019_S11_hi <- data.frame(t(pops.hi))
        colnames(tpop_019_S11_hi) <- c("Seedlings","Saplings","Adults")
        tpop_019_S11_hi$Year <- seq(0,3000,1)   
        tpop_019_S11_lo <- data.frame(t(pops.lo))
        colnames(tpop_019_S11_lo) <- c("Seedlings","Saplings","Adults")
        tpop_019_S11_lo$Year <- seq(0,3000,1)   
        
        pops <- data.frame(N0)       # Part C of the figure: low alpha, perturbed g21
        pops.hi <- data.frame(N0)
        pops.lo <- data.frame(N0)
        for(i in 1:3000){     
          N.old.tmp <- pops[,i]
          N.agg.tmp <- cs*N.old.tmp
          N.old.tmp.hi <- pops.hi[,i]
          N.agg.tmp.hi <- cs*N.old.tmp.hi
          N.old.tmp.lo <- pops.lo[,i]
          N.agg.tmp.lo <- cs*N.old.tmp.lo
          DD.tmp <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp))),0,f13/(1+(alpha.tmp*sum(N.agg.tmp))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          DD.tmp.hi <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.hi))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.hi))),
              boot.ci.g21[2]/(1+(alpha.tmp*sum(N.agg.tmp.hi))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows      
          DD.tmp.lo <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.lo))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.lo))),
              boot.ci.g21[1]/(1+(alpha.tmp*sum(N.agg.tmp.lo))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          N.new.tmp <- DD.tmp%*%N.old.tmp
          N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
          N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
          pops <- data.frame(pops,N.new.tmp)
          pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
          pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
        }
        tpop_019_G21 <- data.frame(t(pops))
        colnames(tpop_019_G21) <- c("Seedlings","Saplings","Adults")
        tpop_019_G21$Year <- seq(0,3000,1)   
        tpop_019_G21_hi <- data.frame(t(pops.hi))
        colnames(tpop_019_G21_hi) <- c("Seedlings","Saplings","Adults")
        tpop_019_G21_hi$Year <- seq(0,3000,1)   
        tpop_019_G21_lo <- data.frame(t(pops.lo))
        colnames(tpop_019_G21_lo) <- c("Seedlings","Saplings","Adults")
        tpop_019_G21_lo$Year <- seq(0,3000,1)   
        
        
        pops <- data.frame(N0)       # Part D of the figure: low alpha, perturbed S22
        pops.hi <- data.frame(N0)
        pops.lo <- data.frame(N0)
        for(i in 1:3000){     
          N.old.tmp <- pops[,i]
          N.agg.tmp <- cs*N.old.tmp
          N.old.tmp.hi <- pops.hi[,i]
          N.agg.tmp.hi <- cs*N.old.tmp.hi
          N.old.tmp.lo <- pops.lo[,i]
          N.agg.tmp.lo <- cs*N.old.tmp.lo
          DD.tmp <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp))),0,f13/(1+(alpha.tmp*sum(N.agg.tmp))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          DD.tmp.hi <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.hi))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.hi))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.hi))),s22*1.01,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows      
          DD.tmp.lo <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.lo))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.lo))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.lo))),s22*.99,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          N.new.tmp <- DD.tmp%*%N.old.tmp
          N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
          N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
          pops <- data.frame(pops,N.new.tmp)
          pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
          pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
        }
        tpop_019_S22 <- data.frame(t(pops))
        colnames(tpop_019_S22) <- c("Seedlings","Saplings","Adults")
        tpop_019_S22$Year <- seq(0,3000,1)   
        tpop_019_S22_hi <- data.frame(t(pops.hi))
        colnames(tpop_019_S22_hi) <- c("Seedlings","Saplings","Adults")
        tpop_019_S22_hi$Year <- seq(0,3000,1)   
        tpop_019_S22_lo <- data.frame(t(pops.lo))
        colnames(tpop_019_S22_lo) <- c("Seedlings","Saplings","Adults")
        tpop_019_S22_lo$Year <- seq(0,3000,1)   
        
        
        pops <- data.frame(N0)       # Part E of the figure: low alpha, perturbed G32
        pops.hi <- data.frame(N0)
        pops.lo <- data.frame(N0)
        for(i in 1:3000){     
          N.old.tmp <- pops[,i]
          N.agg.tmp <- cs*N.old.tmp
          N.old.tmp.hi <- pops.hi[,i]
          N.agg.tmp.hi <- cs*N.old.tmp.hi
          N.old.tmp.lo <- pops.lo[,i]
          N.agg.tmp.lo <- cs*N.old.tmp.lo
          DD.tmp <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp))),0,f13/(1+(alpha.tmp*sum(N.agg.tmp))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          DD.tmp.hi <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.hi))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.hi))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.hi))),s22,0,
              0,g32*1.3,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows      
          DD.tmp.lo <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.lo))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.lo))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.lo))),s22,0,
              0,g32*.7,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          N.new.tmp <- DD.tmp%*%N.old.tmp
          N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
          N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
          pops <- data.frame(pops,N.new.tmp)
          pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
          pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
        }
        tpop_019_G32 <- data.frame(t(pops))
        colnames(tpop_019_G32) <- c("Seedlings","Saplings","Adults")
        tpop_019_G32$Year <- seq(0,3000,1)   
        tpop_019_G32_hi <- data.frame(t(pops.hi))
        colnames(tpop_019_G32_hi) <- c("Seedlings","Saplings","Adults")
        tpop_019_G32_hi$Year <- seq(0,3000,1)   
        tpop_019_G32_lo <- data.frame(t(pops.lo))
        colnames(tpop_019_G32_lo) <- c("Seedlings","Saplings","Adults")
        tpop_019_G32_lo$Year <- seq(0,3000,1)   
        
        
        pops <- data.frame(N0)       # Part F of the figure: low alpha, perturbed S33
        pops.hi <- data.frame(N0)
        pops.lo <- data.frame(N0)
        for(i in 1:3000){     
          N.old.tmp <- pops[,i]
          N.agg.tmp <- cs*N.old.tmp
          N.old.tmp.hi <- pops.hi[,i]
          N.agg.tmp.hi <- cs*N.old.tmp.hi
          N.old.tmp.lo <- pops.lo[,i]
          N.agg.tmp.lo <- cs*N.old.tmp.lo
          DD.tmp <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp))),0,f13/(1+(alpha.tmp*sum(N.agg.tmp))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp))),s22,0,
              0,g32,s33), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          DD.tmp.hi <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.hi))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.hi))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.hi))),s22,0,
              0,g32,boot.ci.s33[2]), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows      
          DD.tmp.lo <- matrix(
            c(s11/(1+(alpha.tmp*sum(N.agg.tmp.lo))),0,(f13)/(1+(alpha.tmp*sum(N.agg.tmp.lo))),
              g21/(1+(alpha.tmp*sum(N.agg.tmp.lo))),s22,0,
              0,g32,boot.ci.s33[1]), # the data elements WITH ALPHAS
            nrow=3,              # number of rows 
            ncol=3,              # number of columns 
            byrow = TRUE)        # fill matrix by rows
          N.new.tmp <- DD.tmp%*%N.old.tmp
          N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
          N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
          pops <- data.frame(pops,N.new.tmp)
          pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
          pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
        }
        tpop_019_S33 <- data.frame(t(pops))
        colnames(tpop_019_S33) <- c("Seedlings","Saplings","Adults")
        tpop_019_S33$Year <- seq(0,3000,1)   
        tpop_019_S33_hi <- data.frame(t(pops.hi))
        colnames(tpop_019_S33_hi) <- c("Seedlings","Saplings","Adults")
        tpop_019_S33_hi$Year <- seq(0,3000,1)   
        tpop_019_S33_lo <- data.frame(t(pops.lo))
        colnames(tpop_019_S33_lo) <- c("Seedlings","Saplings","Adults")
        tpop_019_S33_lo$Year <- seq(0,3000,1)   
        
        
        
        # Plot the population trajectories for Figure 2.1
        par(mfrow = c(2,3))
        plot(tpop_019_F13$Seedlings~tpop_019_F13$Year,col="white",pch = 19,cex=.5,ylim=c(0,100),yaxs="i",xaxs="i",xlab="",ylab="")
        lines(tpop_019_F13$Seedlings~tpop_019_F13$Year,col="forestgreen",cex=3)
        lines(tpop_019_F13_lo$Seedlings~tpop_019_F13_lo$Year,lty=2,col="forestgreen")
        lines(tpop_019_F13_hi$Seedlings~tpop_019_F13_hi$Year,lty=2,col="forestgreen")
        lines(tpop_019_F13$Saplings~tpop_019_F13$Year,col="blue",cex=3)
        lines(tpop_019_F13_lo$Saplings~tpop_019_F13_lo$Year,lty=2,col="blue")
        lines(tpop_019_F13_hi$Saplings~tpop_019_F13_hi$Year,lty=2,col="blue")
        lines(tpop_019_F13$Adults~tpop_019_F13$Year,col="purple",cex=3)
        lines(tpop_019_F13_lo$Adults~tpop_019_F13_lo$Year,lty=2,col="purple")
        lines(tpop_019_F13_hi$Adults~tpop_019_F13_hi$Year,lty=2,col="purple")
        mtext(side=2, line=2, "Individuals/hectare")
        mtext(side=1, line=2, "Years")
        text(1200,95,expression(bold("a ")*italic("F")["1,3"]*" = 1.5 (8.3-12.7)"),cex=1.2) 
        
        plot(tpop_019_S11$Seedlings~tpop_019_S11$Year,col="white",pch = 19,cex=.5,ylim=c(0,100),yaxs="i",xaxs="i",xlab="",ylab="")
        lines(tpop_019_S11$Seedlings~tpop_019_S11$Year,col="forestgreen",cex=3)
        lines(tpop_019_S11_lo$Seedlings~tpop_019_S11_lo$Year,lty=2,col="forestgreen")
        lines(tpop_019_S11_hi$Seedlings~tpop_019_S11_hi$Year,lty=2,col="forestgreen")
        lines(tpop_019_S11$Saplings~tpop_019_S11$Year,col="blue",cex=3)
        lines(tpop_019_S11_lo$Saplings~tpop_019_S11_lo$Year,lty=2,col="blue")
        lines(tpop_019_S11_hi$Saplings~tpop_019_S11_hi$Year,lty=2,col="blue")
        lines(tpop_019_S11$Adults~tpop_019_S11$Year,col="purple",cex=3)
        lines(tpop_019_S11_lo$Adults~tpop_019_S11_lo$Year,lty=2,col="purple")
        lines(tpop_019_S11_hi$Adults~tpop_019_S11_hi$Year,lty=2,col="purple")
        mtext(side=2, line=2, "Individuals/hectare")
        mtext(side=1, line=2, "Years")
        text(1350,95,expression(bold("b ")*italic("S")["1,1"]*" = 0.30 (0.21-0.39)"),cex = 1.2) 
        text(2500,58,"Adults",col="purple",cex=1.5)
        text(2500,20,"Saplings",col="blue",cex=1.5)
        text(2500,38,"Seedlings",col="forestgreen",cex = 1.5)
        
        plot(tpop_019_G21$Seedlings~tpop_019_G21$Year,col="white",pch = 19,cex=.5,ylim=c(0,100),yaxs="i",xaxs="i",xlab="",ylab="")
        lines(tpop_019_G21$Seedlings~tpop_019_G21$Year,col="forestgreen",cex=3)
        lines(tpop_019_G21_lo$Seedlings~tpop_019_G21_lo$Year,lty=2,col="forestgreen")
        lines(tpop_019_G21_hi$Seedlings~tpop_019_G21_hi$Year,lty=2,col="forestgreen")
        lines(tpop_019_G21$Saplings~tpop_019_G21$Year,col="blue",cex=3)
        lines(tpop_019_G21_lo$Saplings~tpop_019_G21_lo$Year,lty=2,col="blue")
        lines(tpop_019_G21_hi$Saplings~tpop_019_G21_hi$Year,lty=2,col="blue")
        lines(tpop_019_G21$Adults~tpop_019_G21$Year,col="purple",cex=3)
        lines(tpop_019_G21_lo$Adults~tpop_019_G21_lo$Year,lty=2,col="purple")
        lines(tpop_019_G21_hi$Adults~tpop_019_G21_hi$Year,lty=2,col="purple")
        mtext(side=2, line=2, "Individuals/hectare")
        mtext(side=1, line=2, "Years")
        text(1500,95,expression(bold("c ")*italic("G")["2,1"]*" = 0.064 (0.024-0.103)"),cex=1.2) 
        
        
        plot(tpop_019_S22$Seedlings~tpop_019_S22$Year,col="white",pch = 19,cex=.5,ylim=c(0,400),yaxs="i",xaxs="i",xlab="",ylab="")
        lines(tpop_019_S22$Seedlings~tpop_019_S22$Year,col="forestgreen",cex=3)
        lines(tpop_019_S22_lo$Seedlings~tpop_019_S22_lo$Year,lty=2,col="forestgreen")
        lines(tpop_019_S22_hi$Seedlings~tpop_019_S22_hi$Year,lty=2,col="forestgreen")
        lines(tpop_019_S22$Saplings~tpop_019_S22$Year,col="blue",cex=3)
        lines(tpop_019_S22_lo$Saplings~tpop_019_S22_lo$Year,lty=2,col="blue")
        lines(tpop_019_S22_hi$Saplings~tpop_019_S22_hi$Year,lty=2,col="blue")
        lines(tpop_019_S22$Adults~tpop_019_S22$Year,col="purple",cex=3)
        lines(tpop_019_S22_lo$Adults~tpop_019_S22_lo$Year,lty=2,col="purple")
        lines(tpop_019_S22_hi$Adults~tpop_019_S22_hi$Year,lty=2,col="purple")
        mtext(side=2, line=2, "Individuals/hectare")
        mtext(side=1, line=2, "Years")
        text(1500,377,expression(bold("d ")*italic("S")["2,2"]*" = 0.990 (0.981-0.999)"),cex=1.2) 
        
        
        plot(tpop_019_G32$Seedlings~tpop_019_G32$Year,col="white",pch = 19,cex=.5,ylim=c(0,100),yaxs="i",xaxs="i",xlab="",ylab="")
        lines(tpop_019_G32$Seedlings~tpop_019_G32$Year,col="forestgreen",cex=3)
        lines(tpop_019_G32_lo$Seedlings~tpop_019_G32_lo$Year,lty=2,col="forestgreen")
        lines(tpop_019_G32_hi$Seedlings~tpop_019_G32_hi$Year,lty=2,col="forestgreen")
        lines(tpop_019_G32$Saplings~tpop_019_G32$Year,col="blue",cex=3)
        lines(tpop_019_G32_lo$Saplings~tpop_019_G32_lo$Year,lty=2,col="blue")
        lines(tpop_019_G32_hi$Saplings~tpop_019_G32_hi$Year,lty=2,col="blue")
        lines(tpop_019_G32$Adults~tpop_019_G32$Year,col="purple",cex=3)
        lines(tpop_019_G32_lo$Adults~tpop_019_G32_lo$Year,lty=2,col="purple")
        lines(tpop_019_G32_hi$Adults~tpop_019_G32_hi$Year,lty=2,col="purple")
        mtext(side=2, line=2, "Individuals/hectare")
        mtext(side=1, line=2, "Years")
        text(1500,95,expression(bold("e ")*italic("G")["3,2"]*" = 0.056 (0.039-0.072)"),cex=1.2) 
        
        plot(tpop_019_S33$Seedlings~tpop_019_S33$Year,col="white",pch = 19,cex=.5,ylim=c(0,100),yaxs="i",xaxs="i",xlab="",ylab="")
        lines(tpop_019_S33$Seedlings~tpop_019_S33$Year,col="forestgreen",cex=3)
        lines(tpop_019_S33_lo$Seedlings~tpop_019_S33_lo$Year,lty=2,col="forestgreen")
        lines(tpop_019_S33_hi$Seedlings~tpop_019_S33_hi$Year,lty=2,col="forestgreen")
        lines(tpop_019_S33$Saplings~tpop_019_S33$Year,col="blue",cex=3)
        lines(tpop_019_S33_lo$Saplings~tpop_019_S33_lo$Year,lty=2,col="blue")
        lines(tpop_019_S33_hi$Saplings~tpop_019_S33_hi$Year,lty=2,col="blue")
        lines(tpop_019_S33$Adults~tpop_019_S33$Year,col="purple",cex=3)
        lines(tpop_019_S33_lo$Adults~tpop_019_S33_lo$Year,lty=2,col="purple")
        lines(tpop_019_S33_hi$Adults~tpop_019_S33_hi$Year,lty=2,col="purple")
        mtext(side=2, line=2, "Individuals/hectare")
        mtext(side=1, line=2, "Years")
        text(1500,95,expression(bold("f ")*italic("S")["3,3"]*" = 0.972 (0.966 - 0.978)"),cex=1.2) 
        
                  #  dev.off()
                    
                    
        # =========== Plot population densities at a stable stage given various values of alpha =============
                    
            ### Conduct a series of simulations of population densities with more values of alpha for Fig 2.2.
                # This figure will have alpha as the x axis
                    # and census density (saplings + adults) as y axis
                    # with solid line for esitmate and dotted lines for boundaries
                    alphas <- seq(.001,.03,.001) # Use more values of alpha
                    
                    # Part A of figure, perturbing f13
                    censuses <- c()
                      censuses.hi <- c()
                      censuses.lo <- c()
                      for(a in 1:length(alphas)){
                            alph.tmp <- alphas[a] # define alpha
                            pops <- data.frame(N0)
                              pops.hi <- data.frame(N0)
                              pops.lo <- data.frame(N0)
                            for(i in 1:500){
                              N.old.tmp <- pops[,i]
                              N.agg.tmp <- cs*N.old.tmp
                              DD.tmp <- matrix(
                                c(s11/(1+(alph.tmp*sum(N.agg.tmp))),0,f13/(1+(alph.tmp*sum(N.agg.tmp))),
                                  g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                                  0,g32,s33), # the data elements WITH ALPHAS
                                nrow=3,              # number of rows 
                                ncol=3,              # number of columns 
                                byrow = TRUE)        # fill matrix by rows
                              N.new.tmp <- DD.tmp%*%N.old.tmp
                              pops <- data.frame(pops,N.new.tmp)
                                        # Now perturb growth rate of interest upward
                                        N.old.tmp.hi <- pops.hi[,i]
                                        N.agg.tmp.hi <- cs*N.old.tmp.hi
                                        DD.tmp.hi <- matrix(
                                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.hi))),0,boot.ci.f13[2]/(1+(alph.tmp*sum(N.agg.tmp.hi))),
                                            g21/(1+(alph.tmp*sum(N.agg.tmp.hi))),s22,0,
                                            0,g32,s33), # the data elements WITH ALPHAS
                                          nrow=3,              # number of rows 
                                          ncol=3,              # number of columns 
                                          byrow = TRUE)        # fill matrix by rows
                                        N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                                        pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
                                        # And perturb growth rate of interest downward
                                        N.old.tmp.lo <- pops.lo[,i]
                                        N.agg.tmp.lo <- cs*N.old.tmp.lo
                                        DD.tmp.lo <- matrix(
                                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.lo))),0,boot.ci.f13[1]/(1+(alph.tmp*sum(N.agg.tmp.lo))),
                                            g21/(1+(alph.tmp*sum(N.agg.tmp.lo))),s22,0,
                                            0,g32,s33), # the data elements WITH ALPHAS
                                          nrow=3,              # number of rows 
                                          ncol=3,              # number of columns 
                                          byrow = TRUE)        # fill matrix by rows
                                        N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                                        pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
                            }
                            censuses <- c(censuses,(pops[2,500]+pops[3,500])) # Save only the final census densities
                              censuses.hi <- c(censuses.hi,(pops.hi[2,500]+pops.hi[3,500]))
                              censuses.lo <- c(censuses.lo,(pops.lo[2,500]+pops.lo[3,500]))
                      }
                    df.f13 <- data.frame(alphas,censuses,censuses.hi,censuses.lo) # densities after X years of growth 
                    
                    
                    # Part B of figure, perturbing s11
                    censuses <- c()
                    censuses.hi <- c()
                    censuses.lo <- c()
                    for(a in 1:length(alphas)){
                      alph.tmp <- alphas[a] # define alpha
                      pops <- data.frame(N0)
                      pops.hi <- data.frame(N0)
                      pops.lo <- data.frame(N0)
                      for(i in 1:500){
                        N.old.tmp <- pops[,i]
                        N.agg.tmp <- cs*N.old.tmp
                        DD.tmp <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp))),0,f13/(1+(alph.tmp*sum(N.agg.tmp))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp <- DD.tmp%*%N.old.tmp
                        pops <- data.frame(pops,N.new.tmp)
                        # Now perturb growth rate of interest upward
                        N.old.tmp.hi <- pops.hi[,i]
                        N.agg.tmp.hi <- cs*N.old.tmp.hi
                        DD.tmp.hi <- matrix(
                          c(boot.ci.s11[2]/(1+(alph.tmp*sum(N.agg.tmp.hi))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.hi))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp.hi))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                        pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
                        # And perturb growth rate of interest downward
                        N.old.tmp.lo <- pops.lo[,i]
                        N.agg.tmp.lo <- cs*N.old.tmp.lo
                        DD.tmp.lo <- matrix(
                          c(boot.ci.s11[1]/(1+(alph.tmp*sum(N.agg.tmp.lo))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.lo))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp.lo))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                        pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
                      }
                      censuses <- c(censuses,(pops[2,500]+pops[3,500])) # Save only the final census densities
                      censuses.hi <- c(censuses.hi,(pops.hi[2,500]+pops.hi[3,500]))
                      censuses.lo <- c(censuses.lo,(pops.lo[2,500]+pops.lo[3,500]))
                    }
                    df.s11 <- data.frame(alphas,censuses,censuses.hi,censuses.lo) # densities after X years of growth 
                    
                    # Part C of figure, perturbing g21
                    censuses <- c()
                    censuses.hi <- c()
                    censuses.lo <- c()
                    for(a in 1:length(alphas)){
                      alph.tmp <- alphas[a] # define alpha
                      pops <- data.frame(N0)
                      pops.hi <- data.frame(N0)
                      pops.lo <- data.frame(N0)
                      for(i in 1:500){
                        N.old.tmp <- pops[,i]
                        N.agg.tmp <- cs*N.old.tmp
                        DD.tmp <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp))),0,f13/(1+(alph.tmp*sum(N.agg.tmp))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp <- DD.tmp%*%N.old.tmp
                        pops <- data.frame(pops,N.new.tmp)
                        # Now perturb growth rate of interest upward
                        N.old.tmp.hi <- pops.hi[,i]
                        N.agg.tmp.hi <- cs*N.old.tmp.hi
                        DD.tmp.hi <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.hi))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.hi))),
                            boot.ci.g21[2]/(1+(alph.tmp*sum(N.agg.tmp.hi))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                        pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
                        # And perturb growth rate of interest downward
                        N.old.tmp.lo <- pops.lo[,i]
                        N.agg.tmp.lo <- cs*N.old.tmp.lo
                        DD.tmp.lo <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.lo))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.lo))),
                            boot.ci.g21[1]/(1+(alph.tmp*sum(N.agg.tmp.lo))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                        pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
                      }
                      censuses <- c(censuses,(pops[2,500]+pops[3,500])) # Save only the final census densities
                      censuses.hi <- c(censuses.hi,(pops.hi[2,500]+pops.hi[3,500]))
                      censuses.lo <- c(censuses.lo,(pops.lo[2,500]+pops.lo[3,500]))
                    }
                    df.g21 <- data.frame(alphas,censuses,censuses.hi,censuses.lo) # densities after X years of growth 
                    
                    # Part D of figure, perturbing s22
                    censuses <- c()
                    censuses.hi <- c()
                    censuses.lo <- c()
                    for(a in 1:length(alphas)){
                      alph.tmp <- alphas[a] # define alpha
                      pops <- data.frame(N0)
                      pops.hi <- data.frame(N0)
                      pops.lo <- data.frame(N0)
                      for(i in 1:500){
                        N.old.tmp <- pops[,i]
                        N.agg.tmp <- cs*N.old.tmp
                        DD.tmp <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp))),0,f13/(1+(alph.tmp*sum(N.agg.tmp))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp <- DD.tmp%*%N.old.tmp
                        pops <- data.frame(pops,N.new.tmp)
                        # Now perturb growth rate of interest upward
                        N.old.tmp.hi <- pops.hi[,i]
                        N.agg.tmp.hi <- cs*N.old.tmp.hi
                        DD.tmp.hi <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.hi))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.hi))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp.hi))),s22*1.01,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                        pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
                        # And perturb growth rate of interest downward
                        N.old.tmp.lo <- pops.lo[,i]
                        N.agg.tmp.lo <- cs*N.old.tmp.lo
                        DD.tmp.lo <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.lo))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.lo))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp.lo))),s22*.99,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                        pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
                      }
                      censuses <- c(censuses,(pops[2,500]+pops[3,500])) # Save only the final census densities
                      censuses.hi <- c(censuses.hi,(pops.hi[2,500]+pops.hi[3,500]))
                      censuses.lo <- c(censuses.lo,(pops.lo[2,500]+pops.lo[3,500]))
                    }
                    df.s22 <- data.frame(alphas,censuses,censuses.hi,censuses.lo) # densities after X years of growth 
                    
                    # Part E of figure, perturbing g32
                    censuses <- c()
                    censuses.hi <- c()
                    censuses.lo <- c()
                    for(a in 1:length(alphas)){
                      alph.tmp <- alphas[a] # define alpha
                      pops <- data.frame(N0)
                      pops.hi <- data.frame(N0)
                      pops.lo <- data.frame(N0)
                      for(i in 1:500){
                        N.old.tmp <- pops[,i]
                        N.agg.tmp <- cs*N.old.tmp
                        DD.tmp <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp))),0,f13/(1+(alph.tmp*sum(N.agg.tmp))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp <- DD.tmp%*%N.old.tmp
                        pops <- data.frame(pops,N.new.tmp)
                        # Now perturb growth rate of interest upward
                        N.old.tmp.hi <- pops.hi[,i]
                        N.agg.tmp.hi <- cs*N.old.tmp.hi
                        DD.tmp.hi <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.hi))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.hi))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp.hi))),s22,0,
                            0,g32*1.3,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                        pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
                        # And perturb growth rate of interest downward
                        N.old.tmp.lo <- pops.lo[,i]
                        N.agg.tmp.lo <- cs*N.old.tmp.lo
                        DD.tmp.lo <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.lo))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.lo))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp.lo))),s22,0,
                            0,g32*0.7,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                        pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
                      }
                      censuses <- c(censuses,(pops[2,500]+pops[3,500])) # Save only the final census densities
                      censuses.hi <- c(censuses.hi,(pops.hi[2,500]+pops.hi[3,500]))
                      censuses.lo <- c(censuses.lo,(pops.lo[2,500]+pops.lo[3,500]))
                    }
                    df.g32 <- data.frame(alphas,censuses,censuses.hi,censuses.lo) # densities after X years of growth 
                    
                    # Part F of figure, perturbing s33
                    censuses <- c()
                    censuses.hi <- c()
                    censuses.lo <- c()
                    for(a in 1:length(alphas)){
                      alph.tmp <- alphas[a] # define alpha
                      pops <- data.frame(N0)
                      pops.hi <- data.frame(N0)
                      pops.lo <- data.frame(N0)
                      for(i in 1:500){
                        N.old.tmp <- pops[,i]
                        N.agg.tmp <- cs*N.old.tmp
                        DD.tmp <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp))),0,f13/(1+(alph.tmp*sum(N.agg.tmp))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp <- DD.tmp%*%N.old.tmp
                        pops <- data.frame(pops,N.new.tmp)
                        # Now perturb growth rate of interest upward
                        N.old.tmp.hi <- pops.hi[,i]
                        N.agg.tmp.hi <- cs*N.old.tmp.hi
                        DD.tmp.hi <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.hi))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.hi))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp.hi))),s22,0,
                            0,g32,boot.ci.s33[2]), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                        pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
                        # And perturb growth rate of interest downward
                        N.old.tmp.lo <- pops.lo[,i]
                        N.agg.tmp.lo <- cs*N.old.tmp.lo
                        DD.tmp.lo <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp.lo))),0,f13/(1+(alph.tmp*sum(N.agg.tmp.lo))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp.lo))),s22,0,
                            0,g32,boot.ci.s33[1]), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                        pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
                      }
                      censuses <- c(censuses,(pops[2,500]+pops[3,500])) # Save only the final census densities
                      censuses.hi <- c(censuses.hi,(pops.hi[2,500]+pops.hi[3,500]))
                      censuses.lo <- c(censuses.lo,(pops.lo[2,500]+pops.lo[3,500]))
                    }
                    df.s33 <- data.frame(alphas,censuses,censuses.hi,censuses.lo) # densities after X years of growth 
                    
                    
                    # Plot the equilibrium population densities at values of alpha for Figure 2.2
                    plot(df.f13$censuses~df.f13$alphas,col="white",pch = 19,cex=.5,ylim=c(0,500),yaxs="i",xaxs="i",xlab="",ylab="")
                    lines(df.f13$censuses~df.f13$alphas,col="black",cex=3)
                    lines(df.f13$censuses.hi~df.f13$alphas,col="black",cex=3,lty=2)
                    lines(df.f13$censuses.lo~df.f13$alphas,col="black",cex=3,lty=2)
                    abline(h=287,col="blue",lty=3)
                    mtext(side=2, line=2, "Individuals/hectare")
                    mtext(side=1, line=2, expression(alpha))
                    text(.012,470,expression(bold("a ")*italic("F")["1,3"]*" = 1.5 (8.3-12.7)"),cex=1.2) 
                    
                    plot(df.s11$censuses~df.s11$alphas,col="white",pch = 19,cex=.5,ylim=c(0,500),yaxs="i",xaxs="i",xlab="",ylab="")
                    lines(df.s11$censuses~df.s11$alphas,col="black",cex=3)
                    lines(df.s11$censuses.hi~df.s11$alphas,col="black",cex=3,lty=2)
                    lines(df.s11$censuses.lo~df.s11$alphas,col="black",cex=3,lty=2)
                    abline(h=287,col="blue",lty=3)
                    mtext(side=2, line=2, "Individuals/hectare")
                    mtext(side=1, line=2, expression(alpha))
                    text(.013,470,expression(bold("b ")*italic("S")["1,1"]*" = 0.30 (0.21-0.39)"),cex=1.2) 
                    
                    plot(df.g21$censuses~df.g21$alphas,col="white",pch = 19,cex=.5,ylim=c(0,500),yaxs="i",xaxs="i",xlab="",ylab="")
                    lines(df.g21$censuses~df.g21$alphas,col="black",cex=3)
                    lines(df.g21$censuses.hi~df.g21$alphas,col="black",cex=3,lty=2)
                    lines(df.g21$censuses.lo~df.g21$alphas,col="black",cex=3,lty=2)
                    abline(h=287,col="blue",lty=3)
                    mtext(side=2, line=2, "Individuals/hectare")
                    mtext(side=1, line=2, expression(alpha))
                    text(.015,470,expression(bold("c ")*italic("G")["2,1"]*" = 0.064 (0.024-0.103)"),cex=1.2) 
                    
                    plot(df.s22$censuses~df.s22$alphas,col="white",pch = 19,cex=.5,ylim=c(0,500),yaxs="i",xaxs="i",xlab="",ylab="")
                    lines(df.s22$censuses~df.s22$alphas,col="black",cex=3)
                    lines(df.s22$censuses.hi~df.s22$alphas,col="black",cex=3,lty=2)
                    lines(df.s22$censuses.lo~df.s22$alphas,col="black",cex=3,lty=2)
                    abline(h=287,col="blue",lty=3)
                    mtext(side=2, line=2, "Individuals/hectare")
                    mtext(side=1, line=2, expression(alpha))
                    text(.015,470,expression(bold("d ")*italic("S")["2,2"]*" = 0.990 (0.981-0.999)"),cex=1.2) 
                    
                    plot(df.g32$censuses~df.g32$alphas,col="white",pch = 19,cex=.5,ylim=c(0,500),yaxs="i",xaxs="i",xlab="",ylab="")
                    lines(df.g32$censuses~df.g32$alphas,col="black",cex=3)
                    lines(df.g32$censuses.hi~df.g32$alphas,col="black",cex=3,lty=2)
                    lines(df.g32$censuses.lo~df.g32$alphas,col="black",cex=3,lty=2)
                    abline(h=287,col="blue",lty=3)
                    mtext(side=2, line=2, "Individuals/hectare")
                    mtext(side=1, line=2, expression(alpha))
                    text(.015,470,expression(bold("e ")*italic("G")["3,2"]*" = 0.056 (0.039-0.072)"),cex=1.2) 
                    
                    plot(df.s33$censuses~df.s33$alphas,col="white",pch = 19,cex=.5,ylim=c(0,500),yaxs="i",xaxs="i",xlab="",ylab="")
                    lines(df.s33$censuses~df.s33$alphas,col="black",cex=3)
                    lines(df.s33$censuses.hi~df.s33$alphas,col="black",cex=3,lty=2)
                    lines(df.s33$censuses.lo~df.s33$alphas,col="black",cex=3,lty=2)
                    abline(h=287,col="blue",lty=3)
                    mtext(side=2, line=2, "Individuals/hectare")
                    mtext(side=1, line=2, expression(alpha))
                    text(.015,470,expression(bold("f ")*italic("S")["3,3"]*" = 0.972 (0.966 - 0.978)"),cex=1.2) 
                    dev.off()
                    # For the actual esimtates, alpha = 0.008 is a useful value
                    
              # Now take a look at how this population would look if all the parameters were 
                    # increased or decreased at the same time
                    
                    alpha.tmp <- .03   # Density dependence factor
                    par(mfrow = c(1,2))
                    pops <- data.frame(N0)       # Part A of the figure: low alpha, perturbed all
                    pops.hi <- data.frame(N0)
                    pops.lo <- data.frame(N0)
                    for(i in 1:3000){     
                      N.old.tmp <- pops[,i]
                      N.agg.tmp <- cs*N.old.tmp
                      N.old.tmp.hi <- pops.hi[,i]
                      N.agg.tmp.hi <- cs*N.old.tmp.hi
                      N.old.tmp.lo <- pops.lo[,i]
                      N.agg.tmp.lo <- cs*N.old.tmp.lo
                      DD.tmp <- matrix(
                        c(s11/(1+(alpha.tmp*sum(N.agg.tmp))),0,f13/(1+(alpha.tmp*sum(N.agg.tmp))),
                          g21/(1+(alpha.tmp*sum(N.agg.tmp))),s22,0,
                          0,g32,s33), # the data elements WITH ALPHAS
                        nrow=3,              # number of rows 
                        ncol=3,              # number of columns 
                        byrow = TRUE)        # fill matrix by rows
                      DD.tmp.hi <- matrix(
                        c(boot.ci.s11[2]/(1+(alpha.tmp*sum(N.agg.tmp.hi))),0,(boot.ci.f13[2])/(1+(alpha.tmp*sum(N.agg.tmp.hi))),
                          boot.ci.g21[2]/(1+(alpha.tmp*sum(N.agg.tmp.hi))),s22*1.01,0,
                          0,g32*1.3,boot.ci.s33[2]), # the data elements WITH ALPHAS
                        nrow=3,              # number of rows 
                        ncol=3,              # number of columns 
                        byrow = TRUE)        # fill matrix by rows      
                      DD.tmp.lo <- matrix(
                        c(boot.ci.s11[1]/(1+(alpha.tmp*sum(N.agg.tmp.hi))),0,(boot.ci.f13[1])/(1+(alpha.tmp*sum(N.agg.tmp.hi))),
                          boot.ci.g21[1]/(1+(alpha.tmp*sum(N.agg.tmp.hi))),s22*.99,0,
                          0,g32*.7,boot.ci.s33[1]), # the data elements WITH ALPHAS
                        nrow=3,              # number of rows 
                        ncol=3,              # number of columns 
                        byrow = TRUE)        # fill matrix by rows
                      N.new.tmp <- DD.tmp%*%N.old.tmp
                      N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                      N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                      pops <- data.frame(pops,N.new.tmp)
                      pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
                      pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
                    }
                    tpop_all <- data.frame(t(pops))
                    colnames(tpop_all) <- c("Seedlings","Saplings","Adults")
                    tpop_all$Year <- seq(0,3000,1)   
                    tpop_all_hi <- data.frame(t(pops.hi))
                    colnames(tpop_all_hi) <- c("Seedlings","Saplings","Adults")
                    tpop_all_hi$Year <- seq(0,3000,1)   
                    tpop_all_lo <- data.frame(t(pops.lo))
                    colnames(tpop_all_lo) <- c("Seedlings","Saplings","Adults")
                    tpop_all_lo$Year <- seq(0,3000,1) 
                    
                    plot(tpop_all$Seedlings~tpop_all$Year,col="white",pch = 19,cex=.5,ylim=c(0,650),yaxs="i",xaxs="i",xlab="",ylab="")
                    lines(tpop_all$Seedlings~tpop_all$Year,col="forestgreen",cex=3)
                    lines(tpop_all_lo$Seedlings~tpop_all_lo$Year,lty=2,col="forestgreen")
                    lines(tpop_all_hi$Seedlings~tpop_all_hi$Year,lty=2,col="forestgreen")
                    lines(tpop_all$Saplings~tpop_all$Year,col="blue",cex=3)
                    lines(tpop_all_lo$Saplings~tpop_all_lo$Year,lty=2,col="blue")
                    lines(tpop_all_hi$Saplings~tpop_all_hi$Year,lty=2,col="blue")
                    lines(tpop_all$Adults~tpop_all$Year,col="purple",cex=3)
                    lines(tpop_all_lo$Adults~tpop_all_lo$Year,lty=2,col="purple")
                    lines(tpop_all_hi$Adults~tpop_all_hi$Year,lty=2,col="purple")
                    mtext(side=2, line=2, "Individuals/hectare")
                    mtext(side=1, line=2, "Years")
                    text(150,625,expression(bold("a "))) 
                    text(2500,525,"Adults",col="purple")
                    text(2500,220,"Saplings",col="blue")
                    text(2500,100,"Seedlings",col="forestgreen")
                    
                    alphas <- seq(0,.1,.001) # Use more values of alpha
                    censuses <- c()
                    censuses.hi <- c()
                    censuses.lo <- c()
                    for(a in 1:length(alphas)){
                      alph.tmp <- alphas[a] # define alpha
                      pops <- data.frame(N0)
                      pops.hi <- data.frame(N0)
                      pops.lo <- data.frame(N0)
                      for(i in 1:3000){
                        N.old.tmp <- pops[,i]
                        N.agg.tmp <- cs*N.old.tmp
                        DD.tmp <- matrix(
                          c(s11/(1+(alph.tmp*sum(N.agg.tmp))),0,f13/(1+(alph.tmp*sum(N.agg.tmp))),
                            g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp <- DD.tmp%*%N.old.tmp
                        pops <- data.frame(pops,N.new.tmp)
                        # Now perturb growth rate of interest upward
                        N.old.tmp.hi <- pops.hi[,i]
                        N.agg.tmp.hi <- cs*N.old.tmp.hi
                        DD.tmp.hi <- matrix(
                          c(boot.ci.s11[2]/(1+(alph.tmp*sum(N.agg.tmp.hi))),0,(boot.ci.f13[2])/(1+(alph.tmp*sum(N.agg.tmp.hi))),
                            boot.ci.g21[2]/(1+(alph.tmp*sum(N.agg.tmp.hi))),s22*1.01,0,
                            0,g32*1.3,boot.ci.s33[2]), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows      
                        N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                        pops.hi <- data.frame(pops.hi,N.new.tmp.hi)
                        # And perturb growth rate of interest downward
                        N.old.tmp.lo <- pops.lo[,i]
                        N.agg.tmp.lo <- cs*N.old.tmp.lo
                        DD.tmp.lo <- matrix(
                          c(boot.ci.s11[1]/(1+(alph.tmp*sum(N.agg.tmp.hi))),0,(boot.ci.f13[1])/(1+(alph.tmp*sum(N.agg.tmp.hi))),
                            boot.ci.g21[1]/(1+(alph.tmp*sum(N.agg.tmp.hi))),s22*.99,0,
                            0,g32*.7,boot.ci.s33[1]), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                        pops.lo <- data.frame(pops.lo,N.new.tmp.lo)
                      }
                      censuses <- c(censuses,(pops[2,300]+pops[3,3000])) # Save only the final census densities
                      censuses.hi <- c(censuses.hi,(pops.hi[2,3000]+pops.hi[3,3000]))
                      censuses.lo <- c(censuses.lo,(pops.lo[2,3000]+pops.lo[3,3000]))
                    }
                    df.all <- data.frame(alphas,censuses,censuses.hi,censuses.lo) # densities after X years of growth 
                    
                    
                    # Plot the equilibrium population densities at values of alpha for Figure 2.3
                    plot(df.all$censuses~df.all$alphas,log="x",col="white",pch = 19,cex=.5,ylim=c(0,600),yaxs="i",xaxs="i",xlab="",ylab="")
                    lines(df.all$censuses~df.all$alphas,col="black",cex=3)
                    lines(df.all$censuses.hi~df.all$alphas,col="black",cex=3,lty=2)
                    lines(df.all$censuses.lo~df.all$alphas,col="black",cex=3,lty=2)
                    abline(h=287,col="blue",lty=3)
                    abline(v=0.0011,col="black",lty=2)
                    mtext(side=2, line=2, "Individuals/hectare")
                    mtext(side=1, line=2, expression(alpha))
                    text(.0012,570,expression(bold("b ")))
                    text(.025,300,"1984-5 census",col="blue")
                    
                    
                    # What's the baseline lambda on the low growth rate matrix?
                    DI.lo <- matrix(
                      c(boot.ci.s11[1],0,(boot.ci.f13[1]),
                        boot.ci.g21[1],(s22*.99),0,
                        0,(g32*.7),boot.ci.s33[1]), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    eigen.analysis(DI.lo)
                    
                    # Baseline lambda on high growth rate matrix:
                    DI.hi <- matrix(
                      c(boot.ci.s11[2],0,(boot.ci.f13[2]),
                        boot.ci.g21[2],(s22*1.01),0,
                        0,(g32*1.3),boot.ci.s33[2]), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    eigen.analysis(DI.hi)
                    
                    
        # ================ Now incorporate the effect of buffel grass on the population: ======================
                  
                    # First, effect of BG on density-independent (DI) pop growth:
                    effect <- rev(seq(.1,1,.1))
                    eigens <- c()
                    eigens.lo <- c()
                    eigens.hi <- c()
                    for(b in 1:length(effect)){
                      DIBG.tmp <- matrix( 
                        c(effect[b]*s11,0,effect[b]*f13,effect[b]*g21,s22,0,0,g32,s33), # the data elements 
                        nrow=3,              # number of rows 
                        ncol=3,              # number of columns 
                        byrow = TRUE)        # fill matrix by rows
                      eigen.tmp <- eigen(DIBG.tmp)[[1]][1] #1.158 
                      eigens <- c(eigens,eigen.tmp)
                            DIBG.lo.tmp <- matrix(
                              c(effect[b]*boot.ci.s11[1],0,(effect[b]*boot.ci.f13[1]),
                                effect[b]*boot.ci.g21[1],(s22*.99),0,
                                0,(g32*.7),boot.ci.s33[1]), # the data elements WITH ALPHAS
                              nrow=3,              # number of rows 
                              ncol=3,              # number of columns 
                              byrow = TRUE)        # fill matrix by rows
                            eigen.lo.tmp <- eigen(DIBG.lo.tmp)[[1]][1] #
                            eigens.lo <- c(eigens.lo,eigen.lo.tmp)
                            DIBG.hi.tmp <- matrix(
                              c(effect[b]*boot.ci.s11[2],0,(effect[b]*boot.ci.f13[2]),
                                effect[b]*boot.ci.g21[2],(s22*1.01),0,
                                0,(g32*1.3),boot.ci.s33[2]), # the data elements WITH ALPHAS
                              nrow=3,              # number of rows 
                              ncol=3,              # number of columns 
                              byrow = TRUE)        # fill matrix by rows
                            eigen.hi.tmp <- eigen(DIBG.hi.tmp)[[1]][1] #
                            eigens.hi <- c(eigens.hi,eigen.hi.tmp)
                    }
                  percent_reduction <- 100*(1-effect)
                  changes_to_di_growth <- data.frame(percent_reduction,eigens,eigens.lo,eigens.hi)
                   # Large effects of BG lead to shrinking populations in the low end of the estimates
                  
                  par(mfrow = c(1,3))
                  #  Plot how eigenvalues change with buffel grass effect for Fig 2.3
                  plot(changes_to_di_growth$eigens~changes_to_di_growth$percent_reduction,col="white",ylim=c(.98,1.35),xlim=c(0,90),xaxs="i",yaxs="i",xlab="",ylab="",cex.axis=1.5)
                  abline(h=1,col="blue",lty=3)
                  lines(changes_to_di_growth$eigens~changes_to_di_growth$percent_reduction,col="black",lwd=3)
                  lines(changes_to_di_growth$eigens.lo~changes_to_di_growth$percent_reduction,col="black",lwd=3,lty=2)
                  lines(changes_to_di_growth$eigens.hi~changes_to_di_growth$percent_reduction,col="black",lwd=3,lty=2)
                  text(30,1.01,"Zero population growth",col="blue",cex=1.5)
                  mtext(side=2, line=2.5, expression(bold(lambda)), cex=1.5)
                  mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                  text(4,1.335,"a",cex=2,font=2)
                  
                  
          # ======================= Calculating effect of buffel grass with density dependence ======================
          # Now plot change in equilibrium densities in the presence of buffel grass for Fig 5b with density dependence
                  alpha <- 0.008
                  alpha.lo <- 0.0004
                  alpha.hi <- 0.107
                  
         # What are euilibrium densities at this alpha in the absence of BG?
                  pops_nobg <- data.frame(N0)
                  pops_nobg.lo <- data.frame(N0)
                  pops_nobg.hi <- data.frame(N0)
                  for(i in 1:5000){
                    N.old.tmp <- pops_nobg[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(s11/(1+(alpha*sum(N.agg.tmp))),0,f13/(1+(alpha*sum(N.agg.tmp))),
                        g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%N.old.tmp
                    pops_nobg <- data.frame(pops_nobg,N.new.tmp)
                    
                    N.old.tmp.lo <- pops_nobg.lo[,i]
                    N.agg.tmp.lo <- cs*N.old.tmp.lo
                    DD.tmp.lo <- matrix(
                      c(boot.ci.s11[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,boot.ci.f13[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                        boot.ci.g21[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22*0.99,0,
                        0,g32*0.7,boot.ci.s33[1]), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                    pops_nobg.lo <- data.frame(pops_nobg.lo,N.new.tmp.lo)
                    
                    N.old.tmp.hi <- pops_nobg.hi[,i]
                    N.agg.tmp.hi <- cs*N.old.tmp.hi
                    DD.tmp.hi <- matrix(
                      c(boot.ci.s11[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,boot.ci.f13[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                        boot.ci.g21[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22*1.01,0,
                        0,g32*1.3,boot.ci.s33[2]), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                    pops_nobg.hi <- data.frame(pops_nobg.hi,N.new.tmp.hi)
                  }
                  # Use these densities to calculate rate of decline for the first centure after BG invasions
                  
                  
                  
             # Calculate equilibrium densities in the presence of buffel grass for Fig 5b-c with density dependence
                    N0b <- pops_nobg[,5001]
                      N0b.lo <- pops_nobg.lo[,5001]
                      N0b.hi <- pops_nobg.hi[,5001]
          
                      eqb_bg <- c()
                      eqb_bg.lo <- c()
                      eqb_bg.hi <- c()
                      bg.century <- c()
                      bg.century.lo <- c()
                      bg.century.hi <- c()
                      bg.effect <- c()
                  for(j in 1:length(effect)){
                    popsb <- data.frame(N0b)
                      popsb.lo <- data.frame(N0b.lo)
                      popsb.hi <- data.frame(N0b.hi)
                    for(i in 1:5000){
                      N.old.tmp <- popsb[,i]
                      N.agg.tmp <- cs*N.old.tmp
                      DD.tmp <- matrix(
                        c(effect[j]*s11/(1+(alpha*sum(N.agg.tmp))),0,effect[j]*f13/(1+(alpha*sum(N.agg.tmp))),
                          effect[j]*g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                          0,g32,s33), # the data elements WITH ALPHAS
                        nrow=3,              # number of rows 
                        ncol=3,              # number of columns 
                        byrow = TRUE)        # fill matrix by rows
                      N.new.tmp <- DD.tmp%*%popsb[,i]
                      popsb <- data.frame(popsb,N.new.tmp)
                      
                      N.old.tmp.lo <- popsb.lo[,i]
                      N.agg.tmp.lo <- cs*N.old.tmp.lo
                      DD.tmp.lo <- matrix(
                        c(effect[j]*boot.ci.s11[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,effect[j]*boot.ci.f13[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                          effect[j]*boot.ci.g21[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22*0.99,0,
                          0,g32*0.7,boot.ci.s33[1]), # the data elements WITH ALPHAS
                        nrow=3,              # number of rows
                        ncol=3,              # number of columns
                        byrow = TRUE)        # fill matrix by rows
                      N.new.tmp.lo <- DD.tmp.lo%*%popsb.lo[,i]
                      popsb.lo <- data.frame(popsb.lo,N.new.tmp.lo)
          
                      N.old.tmp.hi <- popsb.hi[,i]
                      N.agg.tmp.hi <- cs*N.old.tmp.hi
                      DD.tmp.hi <- matrix(
                        c(effect[j]*boot.ci.s11[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,effect[j]*boot.ci.f13[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                          effect[j]*boot.ci.g21[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22*1.01,0,
                          0,g32*1.3,boot.ci.s33[2]), # the data elements WITH ALPHAS
                        nrow=3,              # number of rows
                        ncol=3,              # number of columns
                        byrow = TRUE)        # fill matrix by rows
                      N.new.tmp.hi <- DD.tmp.hi%*%popsb.hi[,i]
                      popsb.hi <- data.frame(popsb.hi,N.new.tmp.hi)
                     }
                   eqb_bg <- c(eqb_bg,(popsb[2,5000]+popsb[3,5000]))
                   eqb_bg.lo <- c(eqb_bg.lo,(popsb.lo[2,5000]+popsb.lo[3,5000]))
                   eqb_bg.hi <- c(eqb_bg.hi,(popsb.hi[2,5000]+popsb.hi[3,5000]))
                   
                   bg.century <- c(bg.century,t(popsb[2,1:100]+popsb[3,1:100]))
                   bg.effect <- c(bg.effect,100*(1-rep(effect[j],100)))
                   bg.century.lo <- c(bg.century.lo,t(popsb.lo[2,1:100]+popsb.lo[3,1:100]))
                   bg.century.hi <- c(bg.century.hi,t(popsb.hi[2,1:100]+popsb.hi[3,1:100]))
                   }
                    # collect densities in each age group
                        
                    eqb_densities_bg <- data.frame(eqb_bg,eqb_bg.hi,eqb_bg.lo,percent_reduction)
                    colnames(eqb_densities_bg) <- c("census","census.hi","census.lo","percent_reduction")
                    
                    # convert them to percentages, to see that each group is reduced by the same percentage
                    eqb_densities_bg$per_red_census <- 100 - round(100*(eqb_densities_bg$census[1]-eqb_densities_bg$census)/(eqb_densities_bg$census[1]),2)
                    eqb_densities_bg$per_red_lo <- 100 - round(100*(eqb_densities_bg$census.lo[1]-eqb_densities_bg$census.lo)/(eqb_densities_bg$census.lo[1]),2)
                    eqb_densities_bg$per_red_hi <- 100 - round(100*(eqb_densities_bg$census.hi[1]-eqb_densities_bg$census.hi)/(eqb_densities_bg$census.hi[1]),2)
                    
                    # Plot population density as % of baseline against the % reduction in recruitment rates (Fig 5b):
                    plot(eqb_densities_bg$per_red_census~eqb_densities_bg$percent_reduction,col="white",xlim=c(0,90),ylim=c(0,108),xaxs="i",yaxs="i",xlab="",ylab="",cex.axis=1.5)
                    lines(eqb_densities_bg$per_red_census~eqb_densities_bg$percent_reduction,col="black",lwd=3)
                    lines(eqb_densities_bg$per_red_lo~eqb_densities_bg$percent_reduction,col="black",lwd=3,lty=2)
                    lines(eqb_densities_bg$per_red_hi~eqb_densities_bg$percent_reduction,col="black",lwd=3,lty=2)
                    mtext(side=2, line=2.5, "% of baseline population", cex=1.1, font = 2)
                    mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                    text(5,104,"b",cex=2, font=2)
                    # 
                      
                    bg.century.per <- 100-(100*((N0b[2]+N0b[3])-bg.century)/(N0b[2]+N0b[3]))
                    bg.century.lo.per <- 100-(100*((N0b.lo[2]+N0b.lo[3])-bg.century.lo)/(N0b.lo[2]+N0b.lo[3]))
                    bg.century.hi.per <- 100-(100*((N0b.hi[2]+N0b.hi[3])-bg.century.hi)/(N0b.hi[2]+N0b.hi[3]))
                    
                    bg.century.ten <- bg.century.per[c(101:200)]
                    bg.century.lo.ten <- bg.century.lo.per[c(101:200)]
                    bg.century.hi.ten <- bg.century.hi.per[c(101:200)]
                    
                    bg.century.thirty <- bg.century.per[c(301:400)]
                    bg.century.lo.thirty <- bg.century.lo.per[c(301:400)]
                    bg.century.hi.thirty <- bg.century.hi.per[c(301:400)]
                    
                    bg.century.ninety <- bg.century.per[c(901:1000)]
                    bg.century.lo.ninety <- bg.century.lo.per[c(901:1000)]
                    bg.century.hi.ninety <- bg.century.hi.per[c(901:1000)]
                    
                    years <- seq(1,100,1)
                    
                    # Plot fig 5c
                    plot(bg.century.ten~years,col="white",xlim=c(0,105),ylim=c(0,108),xaxs="i",yaxs="i",xlab="",ylab="",cex.axis=1.5)
                    lines(bg.century.ten~years,col="black",lwd=4)
                    lines(bg.century.lo.ten~years,col="black",lwd=3,lty=2)
                    lines(bg.century.hi.ten~years,col="black",lwd=3,lty=2)
                    
                    lines(bg.century.thirty~years,col="purple",lwd=4)
                    lines(bg.century.lo.thirty~years,col="purple",lwd=3,lty=2)
                    lines(bg.century.hi.thirty~years,col="purple",lwd=3,lty=2)
                    
                    lines(bg.century.ninety~years,col="tomato",lwd=4)
                    lines(bg.century.lo.ninety~years,col="tomato",lwd=3,lty=2)
                    lines(bg.century.hi.ninety~years,col="tomato",lwd=3,lty=2)
                    mtext(side=2, line=2.5, "% of baseline population", cex=1.1, font = 2)
                    mtext(side=1, line=2.5,"Years", cex=1.1,font = 2)
                    text(5,104,"c",cex=2, font=2)
                    text(96,97,"10%",cex=2, col="black")
                    text(96,70,"30%",cex=2, col="purple")
                    text(96,43,"90%",cex=2, col="tomato")
                    
                  
    # =============== Separate out the modeled degree of variation in each rate separately ===========
                
          # Variability in F13
                # First, effect of BG on density-independent (DI) pop growth:
                effect <- rev(seq(.1,1,.1))
                eigens <- c()
                eigens.lo <- c()
                eigens.hi <- c()
                for(b in 1:length(effect)){
                  DIBG.tmp <- matrix( 
                    c(effect[b]*s11,0,effect[b]*f13,effect[b]*g21,s22,0,0,g32,s33), # the data elements 
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.tmp <- eigen(DIBG.tmp)[[1]][1] #1.158 
                  eigens <- c(eigens,eigen.tmp)
                  DIBG.lo.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*boot.ci.f13[1]),
                      effect[b]*g21,s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.lo.tmp <- eigen(DIBG.lo.tmp)[[1]][1] #
                  eigens.lo <- c(eigens.lo,eigen.lo.tmp)
                  DIBG.hi.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*boot.ci.f13[2]),
                      effect[b]*g21,s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.hi.tmp <- eigen(DIBG.hi.tmp)[[1]][1] #
                  eigens.hi <- c(eigens.hi,eigen.hi.tmp)
                }
                percent_reduction <- 100*(1-effect)
                changes_to_di_growth_f13 <- data.frame(percent_reduction,eigens,eigens.lo,eigens.hi)
                # Large effects of BG lead to shrinking populations in the low end of the estimates
                
                par(mfrow = c(2,3))
                #  Plot how eigenvalues change with buffel grass effect for Fig 2.3
                plot(changes_to_di_growth_f13$eigens~changes_to_di_growth_f13$percent_reduction,col="white",ylim=c(.98,1.35),xlim=c(0,90),xlab="",ylab="",cex.axis=1.5)
                abline(h=1,col="blue",lty=3)
                lines(changes_to_di_growth_f13$eigens~changes_to_di_growth_f13$percent_reduction,col="black",lwd=3)
                lines(changes_to_di_growth_f13$eigens.lo~changes_to_di_growth_f13$percent_reduction,col="black",lwd=3,lty=2)
                lines(changes_to_di_growth_f13$eigens.hi~changes_to_di_growth_f13$percent_reduction,col="black",lwd=3,lty=2)
                text(30,1.01,"Zero population growth",col="blue",cex=1.5)
                mtext(side=2, line=2.5, expression(bold(lambda)), cex=1.5)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                text(29,1.34,expression(bold("a ")*italic("F")["1,3"]*" = 1.5 (8.3-12.7)"),cex=1.5,font=2)
                
               
        # Variability in S11
                # First, effect of BG on density-independent (DI) pop growth:
                effect <- rev(seq(.1,1,.1))
                eigens <- c()
                eigens.lo <- c()
                eigens.hi <- c()
                for(b in 1:length(effect)){
                  DIBG.tmp <- matrix( 
                    c(effect[b]*s11,0,effect[b]*f13,effect[b]*g21,s22,0,0,g32,s33), # the data elements 
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.tmp <- eigen(DIBG.tmp)[[1]][1] #1.158 
                  eigens <- c(eigens,eigen.tmp)
                  DIBG.lo.tmp <- matrix(
                    c(effect[b]*boot.ci.s11[1],0,(effect[b]*f13),
                      effect[b]*g21,s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.lo.tmp <- eigen(DIBG.lo.tmp)[[1]][1] #
                  eigens.lo <- c(eigens.lo,eigen.lo.tmp)
                  DIBG.hi.tmp <- matrix(
                    c(effect[b]*boot.ci.s11[2],0,(effect[b]*f13),
                      effect[b]*g21,s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.hi.tmp <- eigen(DIBG.hi.tmp)[[1]][1] #
                  eigens.hi <- c(eigens.hi,eigen.hi.tmp)
                }
                changes_to_di_growth_s11 <- data.frame(percent_reduction,eigens,eigens.lo,eigens.hi)
                # Large effects of BG lead to shrinking populations in the low end of the estimates
                
                #  Plot how eigenvalues change with buffel grass effect for Fig 2.3
                plot(changes_to_di_growth_s11$eigens~changes_to_di_growth_s11$percent_reduction,col="white",ylim=c(.98,1.35),xlim=c(0,90),xlab="",ylab="",cex.axis=1.5)
                abline(h=1,col="blue",lty=3)
                lines(changes_to_di_growth_s11$eigens~changes_to_di_growth_s11$percent_reduction,col="black",lwd=3)
                lines(changes_to_di_growth_s11$eigens.lo~changes_to_di_growth_s11$percent_reduction,col="black",lwd=3,lty=2)
                lines(changes_to_di_growth_s11$eigens.hi~changes_to_di_growth_s11$percent_reduction,col="black",lwd=3,lty=2)
                text(30,1.01,"Zero population growth",col="blue",cex=1.5)
                mtext(side=2, line=2.5, expression(bold(lambda)), cex=1.5)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                text(33,1.34,expression(bold("b ")*italic("S")["1,1"]*" = 0.30 (0.21-0.39)"),cex=1.5,font=2)
                
      # Variability in G21
                # First, effect of BG on density-independent (DI) pop growth:
                eigens <- c()
                eigens.lo <- c()
                eigens.hi <- c()
                for(b in 1:length(effect)){
                  DIBG.tmp <- matrix( 
                    c(effect[b]*s11,0,effect[b]*f13,effect[b]*g21,s22,0,0,g32,s33), # the data elements 
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.tmp <- eigen(DIBG.tmp)[[1]][1] #1.158 
                  eigens <- c(eigens,eigen.tmp)
                  DIBG.lo.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*f13),
                      effect[b]*boot.ci.g21[1],s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.lo.tmp <- eigen(DIBG.lo.tmp)[[1]][1] #
                  eigens.lo <- c(eigens.lo,eigen.lo.tmp)
                  DIBG.hi.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*f13),
                      effect[b]*boot.ci.g21[2],s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.hi.tmp <- eigen(DIBG.hi.tmp)[[1]][1] #
                  eigens.hi <- c(eigens.hi,eigen.hi.tmp)
                }
                changes_to_di_growth_g21 <- data.frame(percent_reduction,eigens,eigens.lo,eigens.hi)
                # Large effects of BG lead to shrinking populations in the low end of the estimates
                
                #  Plot how eigenvalues change with buffel grass effect for Fig 2.3
                plot(changes_to_di_growth_g21$eigens~changes_to_di_growth_g21$percent_reduction,col="white",ylim=c(.98,1.35),xlim=c(0,90),xlab="",ylab="",cex.axis=1.5)
                abline(h=1,col="blue",lty=3)
                lines(changes_to_di_growth_g21$eigens~changes_to_di_growth_g21$percent_reduction,col="black",lwd=3)
                lines(changes_to_di_growth_g21$eigens.lo~changes_to_di_growth_g21$percent_reduction,col="black",lwd=3,lty=2)
                lines(changes_to_di_growth_g21$eigens.hi~changes_to_di_growth_g21$percent_reduction,col="black",lwd=3,lty=2)
                text(30,1.01,"Zero population growth",col="blue",cex=1.5)
                mtext(side=2, line=2.5, expression(bold(lambda)), cex=1.5)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                text(37,1.34,expression(bold("c ")*italic("G")["2,1"]*" = 0.064 (0.024-0.103)"),cex=1.5,font=2)
                
      # Variability in S22
                # First, effect of BG on density-independent (DI) pop growth:
                eigens <- c()
                eigens.lo <- c()
                eigens.hi <- c()
                for(b in 1:length(effect)){
                  DIBG.tmp <- matrix( 
                    c(effect[b]*s11,0,effect[b]*f13,effect[b]*g21,s22,0,0,g32,s33), # the data elements 
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.tmp <- eigen(DIBG.tmp)[[1]][1] #1.158 
                  eigens <- c(eigens,eigen.tmp)
                  DIBG.lo.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*f13),
                      effect[b]*g21,s22*.99,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.lo.tmp <- eigen(DIBG.lo.tmp)[[1]][1] #
                  eigens.lo <- c(eigens.lo,eigen.lo.tmp)
                  DIBG.hi.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*f13),
                      effect[b]*g21,s22*1.01,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.hi.tmp <- eigen(DIBG.hi.tmp)[[1]][1] #
                  eigens.hi <- c(eigens.hi,eigen.hi.tmp)
                }
                changes_to_di_growth_s22 <- data.frame(percent_reduction,eigens,eigens.lo,eigens.hi)
                # Large effects of BG lead to shrinking populations in the low end of the estimates
                
                #  Plot how eigenvalues change with buffel grass effect for Fig 2.3
                plot(changes_to_di_growth_s22$eigens~changes_to_di_growth_s22$percent_reduction,col="white",ylim=c(.98,1.35),xlim=c(0,90),xlab="",ylab="",cex.axis=1.5)
                abline(h=1,col="blue",lty=3)
                lines(changes_to_di_growth_s22$eigens~changes_to_di_growth_s22$percent_reduction,col="black",lwd=3)
                lines(changes_to_di_growth_s22$eigens.lo~changes_to_di_growth_s22$percent_reduction,col="black",lwd=3,lty=2)
                lines(changes_to_di_growth_s22$eigens.hi~changes_to_di_growth_s22$percent_reduction,col="black",lwd=3,lty=2)
                text(30,1.01,"Zero population growth",col="blue",cex=1.5)
                mtext(side=2, line=2.5, expression(bold(lambda)), cex=1.5)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                text(39,1.34,expression(bold("d ")*italic("S")["2,2"]*" = 0.990 (0.981-0.999)"),cex=1.5,font=2)
    
       # Variability in G32
                # First, effect of BG on density-independent (DI) pop growth:
                eigens <- c()
                eigens.lo <- c()
                eigens.hi <- c()
                for(b in 1:length(effect)){
                  DIBG.tmp <- matrix( 
                    c(effect[b]*s11,0,effect[b]*f13,effect[b]*g21,s22,0,0,g32,s33), # the data elements 
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.tmp <- eigen(DIBG.tmp)[[1]][1] #1.158 
                  eigens <- c(eigens,eigen.tmp)
                  DIBG.lo.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*f13),
                      effect[b]*g21,s22,0,
                      0,g32*.7,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.lo.tmp <- eigen(DIBG.lo.tmp)[[1]][1] #
                  eigens.lo <- c(eigens.lo,eigen.lo.tmp)
                  DIBG.hi.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*f13),
                      effect[b]*g21,s22,0,
                      0,g32*1.3,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.hi.tmp <- eigen(DIBG.hi.tmp)[[1]][1] #
                  eigens.hi <- c(eigens.hi,eigen.hi.tmp)
                }
                changes_to_di_growth_g32 <- data.frame(percent_reduction,eigens,eigens.lo,eigens.hi)
                # Large effects of BG lead to shrinking populations in the low end of the estimates
                
                #  Plot how eigenvalues change with buffel grass effect for Fig 2.3
                plot(changes_to_di_growth_g32$eigens~changes_to_di_growth_g32$percent_reduction,col="white",ylim=c(.98,1.35),xlim=c(0,90),xlab="",ylab="",cex.axis=1.5)
                abline(h=1,col="blue",lty=3)
                lines(changes_to_di_growth_g32$eigens~changes_to_di_growth_g32$percent_reduction,col="black",lwd=3)
                lines(changes_to_di_growth_g32$eigens.lo~changes_to_di_growth_g32$percent_reduction,col="black",lwd=3,lty=2)
                lines(changes_to_di_growth_g32$eigens.hi~changes_to_di_growth_g32$percent_reduction,col="black",lwd=3,lty=2)
                text(30,1.01,"Zero population growth",col="blue",cex=1.5)
                mtext(side=2, line=2.5, expression(bold(lambda)), cex=1.5)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                text(39,1.34,expression(bold("e ")*italic("G")["3,2"]*" = 0.056 (0.039-0.072)"),cex=1.5,font=2)
                
      # Variability in S33
                # First, effect of BG on density-independent (DI) pop growth:
                eigens <- c()
                eigens.lo <- c()
                eigens.hi <- c()
                for(b in 1:length(effect)){
                  DIBG.tmp <- matrix( 
                    c(effect[b]*s11,0,effect[b]*f13,effect[b]*g21,s22,0,0,g32,s33), # the data elements 
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.tmp <- eigen(DIBG.tmp)[[1]][1] #1.158 
                  eigens <- c(eigens,eigen.tmp)
                  DIBG.lo.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*f13),
                      effect[b]*g21,s22,0,
                      0,g32,boot.ci.s33[1]), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.lo.tmp <- eigen(DIBG.lo.tmp)[[1]][1] #
                  eigens.lo <- c(eigens.lo,eigen.lo.tmp)
                  DIBG.hi.tmp <- matrix(
                    c(effect[b]*s11,0,(effect[b]*f13),
                      effect[b]*g21,s22,0,
                      0,g32,boot.ci.s33[2]), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  eigen.hi.tmp <- eigen(DIBG.hi.tmp)[[1]][1] #
                  eigens.hi <- c(eigens.hi,eigen.hi.tmp)
                }
                changes_to_di_growth_s33 <- data.frame(percent_reduction,eigens,eigens.lo,eigens.hi)
                # Large effects of BG lead to shrinking populations in the low end of the estimates
                
                #  Plot how eigenvalues change with buffel grass effect for Fig 2.3
                plot(changes_to_di_growth_s33$eigens~changes_to_di_growth_s33$percent_reduction,col="white",ylim=c(.98,1.35),xlim=c(0,90),xlab="",ylab="",cex.axis=1.5)
                abline(h=1,col="blue",lty=3)
                lines(changes_to_di_growth_s33$eigens~changes_to_di_growth_s33$percent_reduction,col="black",lwd=3)
                lines(changes_to_di_growth_s33$eigens.lo~changes_to_di_growth_s33$percent_reduction,col="black",lwd=3,lty=2)
                lines(changes_to_di_growth_s33$eigens.hi~changes_to_di_growth_s33$percent_reduction,col="black",lwd=3,lty=2)
                text(30,1.01,"Zero population growth",col="blue",cex=1.5)
                mtext(side=2, line=2.5, expression(bold(lambda)), cex=1.5)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                text(40,1.34,expression(bold("f ")*italic("S")["3,3"]*" = 0.972 (0.966 - 0.978)"),cex=1.5,font=2)
                
                
    # ======= And each rate varied separately for % reductionin equilibrium density with density dependence ======            
                  # First varying F13
                  # What are euilibrium densities at this alpha in the absence of BG?
                  pops_nobg <- data.frame(N0)
                  pops_nobg.lo <- data.frame(N0)
                  pops_nobg.hi <- data.frame(N0)
                  for(i in 1:5000){
                    N.old.tmp <- pops_nobg[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(s11/(1+(alpha*sum(N.agg.tmp))),0,f13/(1+(alpha*sum(N.agg.tmp))),
                        g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%N.old.tmp
                    pops_nobg <- data.frame(pops_nobg,N.new.tmp)
                    
                    N.old.tmp.lo <- pops_nobg.lo[,i]
                    N.agg.tmp.lo <- cs*N.old.tmp.lo
                    DD.tmp.lo <- matrix(
                      c(s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,boot.ci.f13[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                        g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                    pops_nobg.lo <- data.frame(pops_nobg.lo,N.new.tmp.lo)
                    
                    N.old.tmp.hi <- pops_nobg.hi[,i]
                    N.agg.tmp.hi <- cs*N.old.tmp.hi
                    DD.tmp.hi <- matrix(
                      c(s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,boot.ci.f13[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                        g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                    pops_nobg.hi <- data.frame(pops_nobg.hi,N.new.tmp.hi)
                  }
                  # Use these densities to calculate rate of decline for the first centure after BG invasions
                  
                    # Calculate equilibrium densities in the presence of buffel grass for Fig 5b-c with density dependence
                    N0b <- pops_nobg[,5000]
                    N0b.lo <- pops_nobg.lo[,5000]
                    N0b.hi <- pops_nobg.hi[,5000]
                    
                    eqb_bg <- c()
                    eqb_bg.lo <- c()
                    eqb_bg.hi <- c()
                    bg.effect <- c()
                    for(j in 1:length(effect)){
                      popsb <- data.frame(N0b)
                      popsb.lo <- data.frame(N0b.lo)
                      popsb.hi <- data.frame(N0b.hi)
                      for(i in 1:5000){
                        N.old.tmp <- popsb[,i]
                        N.agg.tmp <- cs*N.old.tmp
                        DD.tmp <- matrix(
                          c(effect[j]*s11/(1+(alpha*sum(N.agg.tmp))),0,effect[j]*f13/(1+(alpha*sum(N.agg.tmp))),
                            effect[j]*g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows 
                          ncol=3,              # number of columns 
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp <- DD.tmp%*%popsb[,i]
                        popsb <- data.frame(popsb,N.new.tmp)
                        
                        N.old.tmp.lo <- popsb.lo[,i]
                        N.agg.tmp.lo <- cs*N.old.tmp.lo
                        DD.tmp.lo <- matrix(
                          c(effect[j]*s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,effect[j]*boot.ci.f13[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                            effect[j]*g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows
                          ncol=3,              # number of columns
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.lo <- DD.tmp.lo%*%popsb.lo[,i]
                        popsb.lo <- data.frame(popsb.lo,N.new.tmp.lo)
                        
                        N.old.tmp.hi <- popsb.hi[,i]
                        N.agg.tmp.hi <- cs*N.old.tmp.hi
                        DD.tmp.hi <- matrix(
                          c(effect[j]*s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,effect[j]*boot.ci.f13[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                            effect[j]*g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                            0,g32,s33), # the data elements WITH ALPHAS
                          nrow=3,              # number of rows
                          ncol=3,              # number of columns
                          byrow = TRUE)        # fill matrix by rows
                        N.new.tmp.hi <- DD.tmp.hi%*%popsb.hi[,i]
                        popsb.hi <- data.frame(popsb.hi,N.new.tmp.hi)
                      }
                      eqb_bg <- c(eqb_bg,(popsb[2,5000]+popsb[3,5000]))
                      eqb_bg.lo <- c(eqb_bg.lo,(popsb.lo[2,5000]+popsb.lo[3,5000]))
                      eqb_bg.hi <- c(eqb_bg.hi,(popsb.hi[2,5000]+popsb.hi[3,5000]))
                      bg.effect <- c(bg.effect,100*(1-rep(effect[j],100)))
                    }
                  # collect densities in each age group
                  
                  eqb_densities_bg_f13 <- data.frame(eqb_bg,eqb_bg.hi,eqb_bg.lo,percent_reduction)
                  colnames(eqb_densities_bg_f13) <- c("census","census.hi","census.lo","percent_reduction")
                  
                  # convert them to percentages, to see that each group is reduced by the same percentage
                  eqb_densities_bg_f13$per_red_census <- 100 - round(100*(eqb_densities_bg_f13$census[1]-eqb_densities_bg_f13$census)/eqb_densities_bg_f13$census[1],2)
                  eqb_densities_bg_f13$per_red_lo <- 100 - round(100*(eqb_densities_bg_f13$census.lo[1]-eqb_densities_bg_f13$census.lo)/eqb_densities_bg_f13$census.lo[1],2)
                  eqb_densities_bg_f13$per_red_hi <- 100 - round(100*(eqb_densities_bg_f13$census.hi[1]-eqb_densities_bg_f13$census.hi)/eqb_densities_bg_f13$census.hi[1],2)
                  
                  # Plot population density as % of baseline against the % reduction in recruitment rates (Fig 2.4):
                  par(mfrow=c(2,3))
                  plot(eqb_densities_bg_f13$per_red_census~eqb_densities_bg_f13$percent_reduction,col="white",xlim=c(0,90),ylim=c(0,108),xlab="",ylab="",cex.axis=1.5)
                  lines(eqb_densities_bg_f13$per_red_census~eqb_densities_bg_f13$percent_reduction,col="black",lwd=3)
                  lines(eqb_densities_bg_f13$per_red_lo~eqb_densities_bg_f13$percent_reduction,col="black",lwd=3,lty=2)
                  lines(eqb_densities_bg_f13$per_red_hi~eqb_densities_bg_f13$percent_reduction,col="black",lwd=3,lty=2)
                  mtext(side=2, line=2.5, "% of baseline population", cex=1.1, font = 2)
                  mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                  #text(33,106,expression(bold("a ")*italic("F")["1,3"]*" = 1.5 (8.3-12.7)"),cex=2, font=2)
                  text(27,106,expression(bold("a ")*"Seedling emergence"),cex=1.5, font=2)
                  
                # Then varying S11
                # What are euilibrium densities at this alpha in the absence of BG?
                pops_nobg <- data.frame(N0)
                pops_nobg.lo <- data.frame(N0)
                pops_nobg.hi <- data.frame(N0)
                for(i in 1:5000){
                  N.old.tmp <- pops_nobg[,i]
                  N.agg.tmp <- cs*N.old.tmp
                  DD.tmp <- matrix(
                    c(s11/(1+(alpha*sum(N.agg.tmp))),0,f13/(1+(alpha*sum(N.agg.tmp))),
                      g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp <- DD.tmp%*%N.old.tmp
                  pops_nobg <- data.frame(pops_nobg,N.new.tmp)
                  
                  N.old.tmp.lo <- pops_nobg.lo[,i]
                  N.agg.tmp.lo <- cs*N.old.tmp.lo
                  DD.tmp.lo <- matrix(
                    c(boot.ci.s11[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                      g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                  pops_nobg.lo <- data.frame(pops_nobg.lo,N.new.tmp.lo)
                  
                  N.old.tmp.hi <- pops_nobg.hi[,i]
                  N.agg.tmp.hi <- cs*N.old.tmp.hi
                  DD.tmp.hi <- matrix(
                    c(boot.ci.s11[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                      g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                  pops_nobg.hi <- data.frame(pops_nobg.hi,N.new.tmp.hi)
                }
                # Use these densities to calculate rate of decline for the first centure after BG invasions
                
                # Calculate equilibrium densities in the presence of buffel grass for Fig 5b-c with density dependence
                N0b <- pops_nobg[,5001]
                N0b.lo <- pops_nobg.lo[,5001]
                N0b.hi <- pops_nobg.hi[,5001]
                
                eqb_bg <- c()
                eqb_bg.lo <- c()
                eqb_bg.hi <- c()
                bg.effect <- c()
                for(j in 1:length(effect)){
                  popsb <- data.frame(N0b)
                  popsb.lo <- data.frame(N0b.lo)
                  popsb.hi <- data.frame(N0b.hi)
                  for(i in 1:5000){
                    N.old.tmp <- popsb[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(effect[j]*s11/(1+(alpha*sum(N.agg.tmp))),0,effect[j]*f13/(1+(alpha*sum(N.agg.tmp))),
                        effect[j]*g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%popsb[,i]
                    popsb <- data.frame(popsb,N.new.tmp)
                    
                    N.old.tmp.lo <- popsb.lo[,i]
                    N.agg.tmp.lo <- cs*N.old.tmp.lo
                    DD.tmp.lo <- matrix(
                      c(effect[j]*boot.ci.s11[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,effect[j]*f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                        effect[j]*g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.lo <- DD.tmp.lo%*%popsb.lo[,i]
                    popsb.lo <- data.frame(popsb.lo,N.new.tmp.lo)
                    
                    N.old.tmp.hi <- popsb.hi[,i]
                    N.agg.tmp.hi <- cs*N.old.tmp.hi
                    DD.tmp.hi <- matrix(
                      c(effect[j]*boot.ci.s11[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,effect[j]*f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                        effect[j]*g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.hi <- DD.tmp.hi%*%popsb.hi[,i]
                    popsb.hi <- data.frame(popsb.hi,N.new.tmp.hi)
                  }
                  eqb_bg <- c(eqb_bg,(popsb[2,5000]+popsb[3,5000]))
                  eqb_bg.lo <- c(eqb_bg.lo,(popsb.lo[2,5000]+popsb.lo[3,5000]))
                  eqb_bg.hi <- c(eqb_bg.hi,(popsb.hi[2,5000]+popsb.hi[3,5000]))
                  bg.effect <- c(bg.effect,100*(1-rep(effect[j],100)))
                }
                # collect densities in each age group
                
                eqb_densities_bg_s11 <- data.frame(eqb_bg,eqb_bg.hi,eqb_bg.lo,percent_reduction)
                colnames(eqb_densities_bg_s11) <- c("census","census.hi","census.lo","percent_reduction")
                
                # convert them to percentages, to see that each group is reduced by the same percentage
                eqb_densities_bg_s11$per_red_census <- 100 - round(100*(eqb_densities_bg_s11$census[1]-eqb_densities_bg_s11$census)/eqb_densities_bg_s11$census[1],2)
                eqb_densities_bg_s11$per_red_lo <- 100 - round(100*(eqb_densities_bg_s11$census.lo[1]-eqb_densities_bg_s11$census.lo)/eqb_densities_bg_s11$census.lo[1],2)
                eqb_densities_bg_s11$per_red_hi <- 100 - round(100*(eqb_densities_bg_s11$census.hi[1]-eqb_densities_bg_s11$census.hi)/eqb_densities_bg_s11$census.hi[1],2)
                
                # Plot population density as % of baseline against the % reduction in recruitment rates (Fig 2.4):
                plot(eqb_densities_bg_s11$per_red_census~eqb_densities_bg_s11$percent_reduction,col="white",xlim=c(0,90),ylim=c(0,108),xlab="",ylab="",cex.axis=1.5)
                lines(eqb_densities_bg_s11$per_red_census~eqb_densities_bg_s11$percent_reduction,col="black",lwd=3)
                lines(eqb_densities_bg_s11$per_red_lo~eqb_densities_bg_s11$percent_reduction,col="black",lwd=3,lty=2)
                lines(eqb_densities_bg_s11$per_red_hi~eqb_densities_bg_s11$percent_reduction,col="black",lwd=3,lty=2)
                mtext(side=2, line=2.5, "% of baseline population", cex=1.1, font = 2)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                #text(39,106,expression(bold("b ")*italic("S")["1,1"]*" = 0.30 (0.21-0.39)"),cex=2, font=2)       
                text(23,106,expression(bold("b ")*"Seedling survival"),cex=1.5, font=2)
                
                
                # Then varying G21
                # What are euilibrium densities at this alpha in the absence of BG?
                pops_nobg <- data.frame(N0)
                pops_nobg.lo <- data.frame(N0)
                pops_nobg.hi <- data.frame(N0)
                for(i in 1:5000){
                  N.old.tmp <- pops_nobg[,i]
                  N.agg.tmp <- cs*N.old.tmp
                  DD.tmp <- matrix(
                    c(s11/(1+(alpha*sum(N.agg.tmp))),0,f13/(1+(alpha*sum(N.agg.tmp))),
                      g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp <- DD.tmp%*%N.old.tmp
                  pops_nobg <- data.frame(pops_nobg,N.new.tmp)
                  
                  N.old.tmp.lo <- pops_nobg.lo[,i]
                  N.agg.tmp.lo <- cs*N.old.tmp.lo
                  DD.tmp.lo <- matrix(
                    c(s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                      boot.ci.g21[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                  pops_nobg.lo <- data.frame(pops_nobg.lo,N.new.tmp.lo)
                  
                  N.old.tmp.hi <- pops_nobg.hi[,i]
                  N.agg.tmp.hi <- cs*N.old.tmp.hi
                  DD.tmp.hi <- matrix(
                    c(s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                      boot.ci.g21[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                  pops_nobg.hi <- data.frame(pops_nobg.hi,N.new.tmp.hi)
                }
                # Use these densities to calculate rate of decline for the first centure after BG invasions
                
                # Calculate equilibrium densities in the presence of buffel grass for Fig 5b-c with density dependence
                N0b <- pops_nobg[,5001]
                N0b.lo <- pops_nobg.lo[,5001]
                N0b.hi <- pops_nobg.hi[,5001]
                
                eqb_bg <- c()
                eqb_bg.lo <- c()
                eqb_bg.hi <- c()
                bg.effect <- c()
                for(j in 1:length(effect)){
                  popsb <- data.frame(N0b)
                  popsb.lo <- data.frame(N0b.lo)
                  popsb.hi <- data.frame(N0b.hi)
                  for(i in 1:5000){
                    N.old.tmp <- popsb[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(effect[j]*s11/(1+(alpha*sum(N.agg.tmp))),0,effect[j]*f13/(1+(alpha*sum(N.agg.tmp))),
                        effect[j]*g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%popsb[,i]
                    popsb <- data.frame(popsb,N.new.tmp)
                    
                    N.old.tmp.lo <- popsb.lo[,i]
                    N.agg.tmp.lo <- cs*N.old.tmp.lo
                    DD.tmp.lo <- matrix(
                      c(effect[j]*s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,effect[j]*f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                        effect[j]*boot.ci.g21[1]/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.lo <- DD.tmp.lo%*%popsb.lo[,i]
                    popsb.lo <- data.frame(popsb.lo,N.new.tmp.lo)
                    
                    N.old.tmp.hi <- popsb.hi[,i]
                    N.agg.tmp.hi <- cs*N.old.tmp.hi
                    DD.tmp.hi <- matrix(
                      c(effect[j]*s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,effect[j]*f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                        effect[j]*boot.ci.g21[2]/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.hi <- DD.tmp.hi%*%popsb.hi[,i]
                    popsb.hi <- data.frame(popsb.hi,N.new.tmp.hi)
                  }
                  eqb_bg <- c(eqb_bg,(popsb[2,5000]+popsb[3,5000]))
                  eqb_bg.lo <- c(eqb_bg.lo,(popsb.lo[2,5000]+popsb.lo[3,5000]))
                  eqb_bg.hi <- c(eqb_bg.hi,(popsb.hi[2,5000]+popsb.hi[3,5000]))
                  bg.effect <- c(bg.effect,100*(1-rep(effect[j],100)))
                }
                # collect densities in each age group
                
                eqb_densities_bg_g21 <- data.frame(eqb_bg,eqb_bg.hi,eqb_bg.lo,percent_reduction)
                colnames(eqb_densities_bg_g21) <- c("census","census.hi","census.lo","percent_reduction")
                
                # convert them to percentages, to see that each group is reduced by the same percentage
                eqb_densities_bg_g21$per_red_census <- 100 - round(100*(eqb_densities_bg_g21$census[1]-eqb_densities_bg_g21$census)/eqb_densities_bg_g21$census[1],2)
                eqb_densities_bg_g21$per_red_lo <- 100 - round(100*(eqb_densities_bg_g21$census.lo[1]-eqb_densities_bg_g21$census.lo)/eqb_densities_bg_g21$census.lo[1],2)
                eqb_densities_bg_g21$per_red_hi <- 100 - round(100*(eqb_densities_bg_g21$census.hi[1]-eqb_densities_bg_g21$census.hi)/eqb_densities_bg_g21$census.hi[1],2)
                
                # Plot population density as % of baseline against the % reduction in recruitment rates (Fig 2.4):
                plot(eqb_densities_bg_g21$per_red_census~eqb_densities_bg_g21$percent_reduction,col="white",xlim=c(0,90),ylim=c(0,108),xlab="",ylab="",cex.axis=1.5)
                lines(eqb_densities_bg_g21$per_red_census~eqb_densities_bg_g21$percent_reduction,col="black",lwd=3)
                lines(eqb_densities_bg_g21$per_red_lo~eqb_densities_bg_g21$percent_reduction,col="black",lwd=3,lty=2)
                lines(eqb_densities_bg_g21$per_red_hi~eqb_densities_bg_g21$percent_reduction,col="black",lwd=3,lty=2)
                mtext(side=2, line=2.5, "% of baseline population", cex=1.1, font = 2)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
                #text(43,106,expression(bold("c ")*italic("G")["2,1"]*" = 0.064 (0.024-0.103)"),cex=2, font=2)
                text(19,106,expression(bold("c ")*"Establishment"),cex=1.5, font=2)
                
                
                # Then varying S22
                # What are euilibrium densities at this alpha in the absence of BG?
                pops_nobg <- data.frame(N0)
                pops_nobg.lo <- data.frame(N0)
                pops_nobg.hi <- data.frame(N0)
                for(i in 1:5000){
                  N.old.tmp <- pops_nobg[,i]
                  N.agg.tmp <- cs*N.old.tmp
                  DD.tmp <- matrix(
                    c(s11/(1+(alpha*sum(N.agg.tmp))),0,f13/(1+(alpha*sum(N.agg.tmp))),
                      g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp <- DD.tmp%*%N.old.tmp
                  pops_nobg <- data.frame(pops_nobg,N.new.tmp)
                  
                  N.old.tmp.lo <- pops_nobg.lo[,i]
                  N.agg.tmp.lo <- cs*N.old.tmp.lo
                  DD.tmp.lo <- matrix(
                    c(s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                      g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22*.99,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                  pops_nobg.lo <- data.frame(pops_nobg.lo,N.new.tmp.lo)
                  
                  N.old.tmp.hi <- pops_nobg.hi[,i]
                  N.agg.tmp.hi <- cs*N.old.tmp.hi
                  DD.tmp.hi <- matrix(
                    c(s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                      g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22*1.01,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                  pops_nobg.hi <- data.frame(pops_nobg.hi,N.new.tmp.hi)
                }
                # Use these densities to calculate rate of decline for the first centure after BG invasions
                
                # Calculate equilibrium densities in the presence of buffel grass for Fig 5b-c with density dependence
                N0b <- pops_nobg[,5001]
                N0b.lo <- pops_nobg.lo[,5001]
                N0b.hi <- pops_nobg.hi[,5001]
                
                eqb_bg <- c()
                eqb_bg.lo <- c()
                eqb_bg.hi <- c()
                bg.effect <- c()
                for(j in 1:length(effect)){
                  popsb <- data.frame(N0b)
                  popsb.lo <- data.frame(N0b.lo)
                  popsb.hi <- data.frame(N0b.hi)
                  for(i in 1:5000){
                    N.old.tmp <- popsb[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(effect[j]*s11/(1+(alpha*sum(N.agg.tmp))),0,effect[j]*f13/(1+(alpha*sum(N.agg.tmp))),
                        effect[j]*g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%popsb[,i]
                    popsb <- data.frame(popsb,N.new.tmp)
                    
                    N.old.tmp.lo <- popsb.lo[,i]
                    N.agg.tmp.lo <- cs*N.old.tmp.lo
                    DD.tmp.lo <- matrix(
                      c(effect[j]*s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,effect[j]*f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                        effect[j]*g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22*.99,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.lo <- DD.tmp.lo%*%popsb.lo[,i]
                    popsb.lo <- data.frame(popsb.lo,N.new.tmp.lo)
                    
                    N.old.tmp.hi <- popsb.hi[,i]
                    N.agg.tmp.hi <- cs*N.old.tmp.hi
                    DD.tmp.hi <- matrix(
                      c(effect[j]*s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,effect[j]*f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                        effect[j]*g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22*1.01,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.hi <- DD.tmp.hi%*%popsb.hi[,i]
                    popsb.hi <- data.frame(popsb.hi,N.new.tmp.hi)
                  }
                  eqb_bg <- c(eqb_bg,(popsb[2,5000]+popsb[3,5000]))
                  eqb_bg.lo <- c(eqb_bg.lo,(popsb.lo[2,5000]+popsb.lo[3,5000]))
                  eqb_bg.hi <- c(eqb_bg.hi,(popsb.hi[2,5000]+popsb.hi[3,5000]))
                  bg.effect <- c(bg.effect,100*(1-rep(effect[j],100)))
                }
                # collect densities in each age group
                
                eqb_densities_bg_s22 <- data.frame(eqb_bg,eqb_bg.hi,eqb_bg.lo,percent_reduction)
                colnames(eqb_densities_bg_s22) <- c("census","census.hi","census.lo","percent_reduction")
                
                # convert them to percentages, to see that each group is reduced by the same percentage
                eqb_densities_bg_s22$per_red_census <- 100 - round(100*(eqb_densities_bg_s22$census[1]-eqb_densities_bg_s22$census)/eqb_densities_bg_s22$census[1],2)
                eqb_densities_bg_s22$per_red_lo <- 100 - round(100*(eqb_densities_bg_s22$census.lo[1]-eqb_densities_bg_s22$census.lo)/eqb_densities_bg_s22$census.lo[1],2)
                eqb_densities_bg_s22$per_red_hi <- 100 - round(100*(eqb_densities_bg_s22$census.hi[1]-eqb_densities_bg_s22$census.hi)/eqb_densities_bg_s22$census.hi[1],2)
                
                # Plot population density as % of baseline against the % reduction in recruitment rates (Fig 2.4):
                plot(eqb_densities_bg_s22$per_red_census~eqb_densities_bg_s22$percent_reduction,col="white",xlim=c(0,90),ylim=c(0,108),xlab="",ylab="",cex.axis=1.5)
                lines(eqb_densities_bg_s22$per_red_census~eqb_densities_bg_s22$percent_reduction,col="black",lwd=3)
                lines(eqb_densities_bg_s22$per_red_lo~eqb_densities_bg_s22$percent_reduction,col="black",lwd=3,lty=2)
                lines(eqb_densities_bg_s22$per_red_hi~eqb_densities_bg_s22$percent_reduction,col="black",lwd=3,lty=2)
                mtext(side=2, line=2.5, "% of baseline population", cex=1.1, font = 2)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
               # text(45,106,expression(bold("d ")*italic("S")["2,2"]*" = 0.990 (0.981-0.999)"),cex=2, font=2)
                text(23,106,expression(bold("d ")*"Sapling survival"),cex=1.5, font=2)
                
                
                # Then varying G32
                # What are euilibrium densities at this alpha in the absence of BG?
                pops_nobg <- data.frame(N0)
                pops_nobg.lo <- data.frame(N0)
                pops_nobg.hi <- data.frame(N0)
                for(i in 1:5000){
                  N.old.tmp <- pops_nobg[,i]
                  N.agg.tmp <- cs*N.old.tmp
                  DD.tmp <- matrix(
                    c(s11/(1+(alpha*sum(N.agg.tmp))),0,f13/(1+(alpha*sum(N.agg.tmp))),
                      g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp <- DD.tmp%*%N.old.tmp
                  pops_nobg <- data.frame(pops_nobg,N.new.tmp)
                  
                  N.old.tmp.lo <- pops_nobg.lo[,i]
                  N.agg.tmp.lo <- cs*N.old.tmp.lo
                  DD.tmp.lo <- matrix(
                    c(s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                      g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                      0,g32*.7,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                  pops_nobg.lo <- data.frame(pops_nobg.lo,N.new.tmp.lo)
                  
                  N.old.tmp.hi <- pops_nobg.hi[,i]
                  N.agg.tmp.hi <- cs*N.old.tmp.hi
                  DD.tmp.hi <- matrix(
                    c(s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                      g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                      0,g32*1.3,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                  pops_nobg.hi <- data.frame(pops_nobg.hi,N.new.tmp.hi)
                }
                # Use these densities to calculate rate of decline for the first centure after BG invasions
                
                # Calculate equilibrium densities in the presence of buffel grass for Fig 5b-c with density dependence
                N0b <- pops_nobg[,5001]
                N0b.lo <- pops_nobg.lo[,5001]
                N0b.hi <- pops_nobg.hi[,5001]
                
                eqb_bg <- c()
                eqb_bg.lo <- c()
                eqb_bg.hi <- c()
                bg.effect <- c()
                for(j in 1:length(effect)){
                  popsb <- data.frame(N0b)
                  popsb.lo <- data.frame(N0b.lo)
                  popsb.hi <- data.frame(N0b.hi)
                  for(i in 1:5000){
                    N.old.tmp <- popsb[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(effect[j]*s11/(1+(alpha*sum(N.agg.tmp))),0,effect[j]*f13/(1+(alpha*sum(N.agg.tmp))),
                        effect[j]*g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%popsb[,i]
                    popsb <- data.frame(popsb,N.new.tmp)
                    
                    N.old.tmp.lo <- popsb.lo[,i]
                    N.agg.tmp.lo <- cs*N.old.tmp.lo
                    DD.tmp.lo <- matrix(
                      c(effect[j]*s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,effect[j]*f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                        effect[j]*g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                        0,g32*.7,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.lo <- DD.tmp.lo%*%popsb.lo[,i]
                    popsb.lo <- data.frame(popsb.lo,N.new.tmp.lo)
                    
                    N.old.tmp.hi <- popsb.hi[,i]
                    N.agg.tmp.hi <- cs*N.old.tmp.hi
                    DD.tmp.hi <- matrix(
                      c(effect[j]*s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,effect[j]*f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                        effect[j]*g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                        0,g32*1.3,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.hi <- DD.tmp.hi%*%popsb.hi[,i]
                    popsb.hi <- data.frame(popsb.hi,N.new.tmp.hi)
                  }
                  eqb_bg <- c(eqb_bg,(popsb[2,5000]+popsb[3,5000]))
                  eqb_bg.lo <- c(eqb_bg.lo,(popsb.lo[2,5000]+popsb.lo[3,5000]))
                  eqb_bg.hi <- c(eqb_bg.hi,(popsb.hi[2,5000]+popsb.hi[3,5000]))
                  bg.effect <- c(bg.effect,100*(1-rep(effect[j],100)))
                }
                # collect densities in each age group
                
                eqb_densities_bg_g32 <- data.frame(eqb_bg,eqb_bg.hi,eqb_bg.lo,percent_reduction)
                colnames(eqb_densities_bg_g32) <- c("census","census.hi","census.lo","percent_reduction")
                
                # convert them to percentages, to see that each group is reduced by the same percentage
                eqb_densities_bg_g32$per_red_census <- 100 - round(100*(eqb_densities_bg_g32$census[1]-eqb_densities_bg_g32$census)/eqb_densities_bg_g32$census[1],2)
                eqb_densities_bg_g32$per_red_lo <- 100 - round(100*(eqb_densities_bg_g32$census.lo[1]-eqb_densities_bg_g32$census.lo)/eqb_densities_bg_g32$census.lo[1],2)
                eqb_densities_bg_g32$per_red_hi <- 100 - round(100*(eqb_densities_bg_g32$census.hi[1]-eqb_densities_bg_g32$census.hi)/eqb_densities_bg_g32$census.hi[1],2)
                
                # Plot population density as % of baseline against the % reduction in recruitment rates (Fig 2.4):
                plot(eqb_densities_bg_g32$per_red_census~eqb_densities_bg_g32$percent_reduction,col="white",xlim=c(0,90),ylim=c(0,108),xlab="",ylab="",cex.axis=1.5)
                lines(eqb_densities_bg_g32$per_red_census~eqb_densities_bg_g32$percent_reduction,col="black",lwd=3)
                lines(eqb_densities_bg_g32$per_red_lo~eqb_densities_bg_g32$percent_reduction,col="black",lwd=3,lty=2)
                lines(eqb_densities_bg_g32$per_red_hi~eqb_densities_bg_g32$percent_reduction,col="black",lwd=3,lty=2)
                mtext(side=2, line=2.5, "% of baseline population", cex=1.1, font = 2)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
               # text(45,106,expression(bold("e ")*italic("G")["3,2"]*" = 0.056 (0.039-0.072)"),cex=2, font=2)
                text(25,106,expression(bold("e ")*"Age of reproduction"),cex=1.5, font=2)
                
                # Then varying S33
                # What are euilibrium densities at this alpha in the absence of BG?
                pops_nobg <- data.frame(N0)
                pops_nobg.lo <- data.frame(N0)
                pops_nobg.hi <- data.frame(N0)
                for(i in 1:5000){
                  N.old.tmp <- pops_nobg[,i]
                  N.agg.tmp <- cs*N.old.tmp
                  DD.tmp <- matrix(
                    c(s11/(1+(alpha*sum(N.agg.tmp))),0,f13/(1+(alpha*sum(N.agg.tmp))),
                      g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                      0,g32,s33), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp <- DD.tmp%*%N.old.tmp
                  pops_nobg <- data.frame(pops_nobg,N.new.tmp)
                  
                  N.old.tmp.lo <- pops_nobg.lo[,i]
                  N.agg.tmp.lo <- cs*N.old.tmp.lo
                  DD.tmp.lo <- matrix(
                    c(s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                      g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                      0,g32,boot.ci.s33[1]), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.lo <- DD.tmp.lo%*%N.old.tmp.lo
                  pops_nobg.lo <- data.frame(pops_nobg.lo,N.new.tmp.lo)
                  
                  N.old.tmp.hi <- pops_nobg.hi[,i]
                  N.agg.tmp.hi <- cs*N.old.tmp.hi
                  DD.tmp.hi <- matrix(
                    c(s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                      g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                      0,g32,boot.ci.s33[2]), # the data elements WITH ALPHAS
                    nrow=3,              # number of rows 
                    ncol=3,              # number of columns 
                    byrow = TRUE)        # fill matrix by rows
                  N.new.tmp.hi <- DD.tmp.hi%*%N.old.tmp.hi
                  pops_nobg.hi <- data.frame(pops_nobg.hi,N.new.tmp.hi)
                }
                # Use these densities to calculate rate of decline for the first centure after BG invasions
                
                # Calculate equilibrium densities in the presence of buffel grass for Fig 5b-c with density dependence
                N0b <- pops_nobg[,5001]
                N0b.lo <- pops_nobg.lo[,5001]
                N0b.hi <- pops_nobg.hi[,5001]
                
                eqb_bg <- c()
                eqb_bg.lo <- c()
                eqb_bg.hi <- c()
                bg.effect <- c()
                for(j in 1:length(effect)){
                  popsb <- data.frame(N0b)
                  popsb.lo <- data.frame(N0b.lo)
                  popsb.hi <- data.frame(N0b.hi)
                  for(i in 1:5000){
                    N.old.tmp <- popsb[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(effect[j]*s11/(1+(alpha*sum(N.agg.tmp))),0,effect[j]*f13/(1+(alpha*sum(N.agg.tmp))),
                        effect[j]*g21/(1+(alpha*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%popsb[,i]
                    popsb <- data.frame(popsb,N.new.tmp)
                    
                    N.old.tmp.lo <- popsb.lo[,i]
                    N.agg.tmp.lo <- cs*N.old.tmp.lo
                    DD.tmp.lo <- matrix(
                      c(effect[j]*s11/(1+(alpha.lo*sum(N.agg.tmp.lo))),0,effect[j]*f13/(1+(alpha.lo*sum(N.agg.tmp.lo))),
                        effect[j]*g21/(1+(alpha.lo*sum(N.agg.tmp.lo))),s22,0,
                        0,g32,boot.ci.s33[1]), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.lo <- DD.tmp.lo%*%popsb.lo[,i]
                    popsb.lo <- data.frame(popsb.lo,N.new.tmp.lo)
                    
                    N.old.tmp.hi <- popsb.hi[,i]
                    N.agg.tmp.hi <- cs*N.old.tmp.hi
                    DD.tmp.hi <- matrix(
                      c(effect[j]*s11/(1+(alpha.hi*sum(N.agg.tmp.hi))),0,effect[j]*f13/(1+(alpha.hi*sum(N.agg.tmp.hi))),
                        effect[j]*g21/(1+(alpha.hi*sum(N.agg.tmp.hi))),s22,0,
                        0,g32,boot.ci.s33[2]), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows
                      ncol=3,              # number of columns
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp.hi <- DD.tmp.hi%*%popsb.hi[,i]
                    popsb.hi <- data.frame(popsb.hi,N.new.tmp.hi)
                  }
                  eqb_bg <- c(eqb_bg,(popsb[2,5000]+popsb[3,5000]))
                  eqb_bg.lo <- c(eqb_bg.lo,(popsb.lo[2,5000]+popsb.lo[3,5000]))
                  eqb_bg.hi <- c(eqb_bg.hi,(popsb.hi[2,5000]+popsb.hi[3,5000]))
                  bg.effect <- c(bg.effect,100*(1-rep(effect[j],100)))
                }
                # collect densities in each age group
                
                eqb_densities_bg_s33 <- data.frame(eqb_bg,eqb_bg.hi,eqb_bg.lo,percent_reduction)
                colnames(eqb_densities_bg_s33) <- c("census","census.hi","census.lo","percent_reduction")
                
                # convert them to percentages, to see that each group is reduced by the same percentage
                eqb_densities_bg_s33$per_red_census <- 100 - round(100*(eqb_densities_bg_s33$census[1]-eqb_densities_bg_s33$census)/eqb_densities_bg_s33$census[1],2)
                eqb_densities_bg_s33$per_red_lo <- 100 - round(100*(eqb_densities_bg_s33$census.lo[1]-eqb_densities_bg_s33$census.lo)/eqb_densities_bg_s33$census.lo[1],2)
                eqb_densities_bg_s33$per_red_hi <- 100 - round(100*(eqb_densities_bg_s33$census.hi[1]-eqb_densities_bg_s33$census.hi)/eqb_densities_bg_s33$census.hi[1],2)
                
                # Plot population density as % of baseline against the % reduction in recruitment rates (Fig 2.4):
                plot(eqb_densities_bg_s33$per_red_census~eqb_densities_bg_s33$percent_reduction,col="white",xlim=c(0,90),ylim=c(0,108),xlab="",ylab="",cex.axis=1.5)
                lines(eqb_densities_bg_s33$per_red_census~eqb_densities_bg_s33$percent_reduction,col="black",lwd=3)
                lines(eqb_densities_bg_s33$per_red_lo~eqb_densities_bg_s33$percent_reduction,col="black",lwd=3,lty=2)
                lines(eqb_densities_bg_s33$per_red_hi~eqb_densities_bg_s33$percent_reduction,col="black",lwd=3,lty=2)
                mtext(side=2, line=2.5, "% of baseline population", cex=1.1, font = 2)
                mtext(side=1, line=2.5,"% reduction in recruitment rates", cex=1.1,font = 2)
               # text(45,106,expression(bold("f ")*italic("S")["3,3"]*" = 0.972 (0.966 - 0.978)"),cex=2, font=2)
                text(17,106,expression(bold("f ")*"Adult survival"),cex=1.5, font=2)
                
                
                
                
                
    # ============= Illustration that alpha does not change % effect of buffel grass ===============
                # Repeating growth simulation with only best estimate of population parameters 
                # for varying values of alpha to demonstrate its value does not matter for 
                # relative effect of buffel grass
                alphas <- seq(0,.1,.01) # Use more values of alpha
                censuses <- c()
                structures <- data.frame(N0)
                for(a in 1:length(alphas)){
                  alph.tmp <- alphas[a] # define alpha
                  pops <- data.frame(N0)
                  for(i in 1:3000){
                    N.old.tmp <- pops[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(s11/(1+(alph.tmp*sum(N.agg.tmp))),0,f13/(1+(alph.tmp*sum(N.agg.tmp))),
                        g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%N.old.tmp
                    pops <- data.frame(pops,N.new.tmp)
                  }
                  censuses <- c(censuses,(pops[2,3000]+pops[3,3000])) # Save only the final census densities
                  structures <- cbind(structures,pops[,3000])
                }
                df.all <- data.frame(alphas,censuses) # densities after 3000 years of growth 
                structures <- structures[,-1]
                eqb_bg <- c()
                for(j in 1:length(alphas)){
                  alph.tmp <- alphas[j] # define alpha
                  struc.tmp <- structures[,j] # eqb population w/o buffel to start from
                  popsb <- data.frame(struc.tmp)
                  for(i in 1:1000){
                    N.old.tmp <- popsb[,i]
                    N.agg.tmp <- cs*N.old.tmp
                    DD.tmp <- matrix(
                      c(.5*s11/(1+(alph.tmp*sum(N.agg.tmp))),0,.5*f13/(1+(alph.tmp*sum(N.agg.tmp))),
                        .5*g21/(1+(alph.tmp*sum(N.agg.tmp))),s22,0,
                        0,g32,s33), # the data elements WITH ALPHAS
                      nrow=3,              # number of rows 
                      ncol=3,              # number of columns 
                      byrow = TRUE)        # fill matrix by rows
                    N.new.tmp <- DD.tmp%*%popsb[,i]
                    popsb <- data.frame(popsb,N.new.tmp)
                  }
                  eqb_bg <- c(eqb_bg,(popsb[2,1000]+popsb[3,1000]))
                }
                # collect densities in each age group
                eqb_densities_bg <- data.frame(eqb_bg,alphas)
                colnames(eqb_densities_bg) <- c("census","alpha")
                eqb_densities_bg$nobg <- df.all$censuses
                eqb_densities_bg <- eqb_densities_bg[-1,]
                # convert them to percentages, to see that each group is reduced by the same percentage
                eqb_densities_bg$percent_reduction <- 100-round(100*(eqb_densities_bg$nobg-eqb_densities_bg$census)/eqb_densities_bg$nobg,2)
    # Alpha value does NOT chnage the percent reduction in palo verde population caused by buffel grass in this model
              