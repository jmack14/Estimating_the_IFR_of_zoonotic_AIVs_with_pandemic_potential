# ============================================================
# Sensitivity analysis: genotype weighting assumptions
# ============================================================

# ------------------------------------------------------------
# Generate genotype samples once
# ------------------------------------------------------------

generate_genotypes <- function(param_ranges, num_samples){
  
  lhs_samples <- randomLHS(num_samples, length(param_ranges))
  
  transformed_samples <- data.frame(
    
    a = qunif(lhs_samples[,1],
              min = param_ranges$a[1],
              max = param_ranges$a[2]),
    
    Rstar = qtriangle(lhs_samples[,2],
                      a = param_ranges$Rstar[1],
                      b = param_ranges$Rstar[2],
                      c = param_ranges$Rstar[3]),
    
    R0 = qexp(lhs_samples[,3],
              rate = param_ranges$R0$rate),
    
    psi = qgamma(lhs_samples[,4],
                 shape = param_ranges$psi$shape,
                 scale = param_ranges$psi$scale)
  )
  
  
  # pandemic probability for each genotype
  transformed_samples$pi <- apply(
    transformed_samples[,c("a","Rstar","R0","psi")],
    1,
    model
  )
  
  transformed_samples
  
}


# ------------------------------------------------------------
# Calculate infections under a weighting scheme
# ------------------------------------------------------------

run_weighting <- function(genotypes, weighting){
  
  n <- nrow(genotypes)
  
  
  # genotype weights rho_i
  if(weighting == "equal"){
    
    rho <- rep(1/n,n)
    
  }
  
  
  if(weighting == "Dirichlet"){
    
    rho <- rgamma(n,shape=1)
    rho <- rho/sum(rho)
    
  }
  
  
  if(weighting == "R0_positive"){
    
    rho <- genotypes$R0
    rho <- rho/sum(rho)
    
  }
  
  
  if(weighting == "R0_negative"){
    
    rho <- 1-genotypes$R0
    rho <- rho/sum(rho)
    
  }
  
  
  # weighted probability of pandemic per infection
  bar_pi <- sum(rho * genotypes$pi)
  
  
  # solve for zoonotic spillovers
  Lambda <- function(mean_nz, Lambda_target){
    
    nz <- 0:10000
    
    PNz_nz <- dpois(nz,mean_nz)
    
    1 - sum(PNz_nz*(1-bar_pi)^nz) - Lambda_target
    
  }
  
  
  mean_nz <- numeric(n)
  
  
  for(i in seq_len(n)){
    
    mean_nz[i] <- uniroot(
      Lambda,
      interval=c(0,1e6),
      Lambda_target=1/genotypes$psi[i]
    )$root
    
  }
  
  
  # weighted mean infections per spillover
  weighted_m <- sum(
    rho *
      (1/(1-genotypes$R0))
  )
  
  
  # annual infections
  mean_nh <- weighted_m * mean_nz
  
  
  # IFR
  IFR <- 20.6 / mean_nh
  
  
  data.frame(
    weighting = weighting,
    mean_nh = mean_nh,
    IFR = IFR
  )
  
}



# ============================================================
# Run sensitivity analysis
# ============================================================

num_samples <- 1000


param_ranges <- list(
  
  R0 = list(rate=20),
  
  a = c(0.000012,0.000024),
  
  Rstar = c(1,2,1.1),
  
  psi = list(
    shape=delta/theta0,
    scale=theta0
  )
  
)



# Generate one common genotype sample
genotypes <- generate_genotypes(
  param_ranges,
  num_samples
)



# Apply alternative genotype weights
SA_results <- bind_rows(
  
  run_weighting(genotypes,"equal"),
  
  run_weighting(genotypes,"Dirichlet"),
  
  run_weighting(genotypes,"R0_positive"),
  
  run_weighting(genotypes,"R0_negative")
  
)



# ============================================================
# Summary table
# ============================================================

SA_summary <- SA_results %>%
  
  group_by(weighting) %>%
  
  summarise(
    
    median_infections =
      median(mean_nh),
    
    infections_lower95 =
      quantile(mean_nh,0.025),
    
    infections_upper95 =
      quantile(mean_nh,0.975),
    
    
    median_IFR =
      20.6 / median(mean_nh),
    
    IFR_lower95 =
      20.6 / quantile(mean_nh,0.975),
    
    IFR_upper95 =
      20.6 / quantile(mean_nh,0.025)
    
  )


print(SA_summary)


write.csv(
  SA_summary,
  "Genotype_weighting_sensitivity_analysis.csv",
  row.names=FALSE
)