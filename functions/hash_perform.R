

###########################################################################
###  Hashimoto Vulnerability Function
###########################################################################
hash_perform <- function(demand_ts, shortage_col){
### Follows Hashimoto, T., Stedinger, J.R., Loucks, D.P., 1982. Reliability, resiliency, and vulnerability criteria for water resource system performance evaluation. Water Resour. Res. 18, 14–20. https://doi.org/10.1029/WR018i001p00014
### Also used https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2002WR001778

shortage <- demand_ts %>% select(shortage_col)

sat <- shortage <= 0
unsat <- shortage > 0

sat_lag <- c(sat[seq(2, length(sat))], NA)
w_trans <- unsat == TRUE & sat_lag == TRUE

time_length <- sum(!is.na(shortage))

### Reliability
### Fraction of time that demand is met
reliability <- sum(sat, na.rm=TRUE)/time_length

### Resilience
### Average probability of recovery in any given failure time step
resilience <- sum(w_trans, na.rm=TRUE)/ (time_length - sum(sat, na.rm=TRUE))

### Vulnerability
### Maximum shortage
vulnerability <- max(shortage, na.rm=TRUE)

return(list(reliability=reliability, resilience=resilience, vulnerability = vulnerability))
}

