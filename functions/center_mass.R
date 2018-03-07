		center_mass <- function (weight, distance) {
			moment <- weight * distance

			moment_sum <- sum(moment)
			weight_sum <- sum(weight)
			return(moment_sum/weight_sum)
		}
		