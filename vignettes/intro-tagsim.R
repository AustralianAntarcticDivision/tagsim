## ------------------------------------------------------------------------
library(tagsim)

## ---- eval=FALSE---------------------------------------------------------
#  ?tagsim

## ------------------------------------------------------------------------
pop <- list("initial" = 1e7, "nat_mort" = 0.1)

## ------------------------------------------------------------------------
rec <- list("type" = "constant", "initial" = 1e6, "mu"=0, "s"=1,
            "spat_dist"="uniform", "stochastic_rec"=TRUE)

## ------------------------------------------------------------------------
harvest <- list("type" = "const_exploit", "rate" = 0.05, "max" = 0.95,
                "ricker"=1)

## ------------------------------------------------------------------------
tag <- list("rel_rate" = 0.05, "mort" = 0, "shed" = 0, "report" = 1)

## ------------------------------------------------------------------------
fish <- list("type" = "random", "prop" = 0.5)

## ------------------------------------------------------------------------
move <- list("type" = "complete", "prob" = 0.05)

## ------------------------------------------------------------------------
assess <- list("type" = "survey", "cv"=0.2)

## ------------------------------------------------------------------------
tmp <- create_control(years = 1:5, 
                      regions = 1:10, 
                      pop_pars = pop, 
                      rec_pars = rec,
                      harvest_pars = harvest, 
                      tag_pars = tag, 
                      fish_pars = fish, 
                      move_pars = move, 
                      assess_pars = assess)

## ------------------------------------------------------------------------
str(tmp)

## ------------------------------------------------------------------------
pop <- list("initial" = 1e7, "nat_mort" = 0.1)
rec <- list("type" = "constant", "initial" = 5e5, "mu"=0, "s"=1,
            "spat_dist"="uniform", "stochastic_rec"=TRUE)
tag <- list("rel_rate" = 0.05, "mort" = 0, "shed" = 0, "report" = 1)
fish <- list("type" = "random", "prop" = 1)
move <- list("type" = "complete", "prob" = 0)
assess <- list("type" = "survey", "cv"=0.1)


## ------------------------------------------------------------------------
reps <- 500
exploit_rate <- 0.025 * 1:10
res <- matrix(0, reps, length(exploit_rate))

## ------------------------------------------------------------------------
for(i in 1:length(exploit_rate)){
  ## update the harvest rate parameters
  harvest <- list("type" = "const_exploit", "rate" = exploit_rate[i], 
                  "max" = 0.95, "ricker"=1)
  ## create the control file
  ctrl <- create_control(years = 1:5, 
                         regions = 1, 
                         pop_pars = pop, 
                         rec_pars = rec,
                         harvest_pars = harvest, 
                         tag_pars = tag, 
                         fish_pars = fish, 
                         move_pars = move, 
                         assess_pars = assess)
  ## run the ith simulation
  obj <- run_sim(ctrl, reps)
  ## save the end of year % unfished biomass
  res[, i] <- obj$true_N[6,] / obj$true_N[1,]
}


## ------------------------------------------------------------------------
mns <- colMeans(res)

## ------------------------------------------------------------------------
plot(mns ~ exploit_rate, las=2, ylim=c(0,1), type="b",
     main="Biomass after 5 years of fishing",
        ylab="Biomass as % unfished", xlab="Exploitation rate")


## ------------------------------------------------------------------------
sessionInfo()

