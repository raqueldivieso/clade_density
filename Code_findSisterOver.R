findSisterOver <- function(tree, overMatrix) {
	library(phytools)
	sisters <- matrix(NA, ncol = 2, nrow = length(tree$tip.label))
	sisters[, 1] <- tree$tip.label
	for (i in 1:length(tree$tip.label)) {
		tmp <- phytools::getSisters(tree, tree$tip.label[i], mode = "label")
		if (!is.null(tmp$tips)) {
			sisters[i, 2] <- tmp$tips
		}
	}
	sisters <- sisters[-which(is.na(sisters[, 2])), ]

	for (i in 1:dim(sisters)[1]) {
		sisters[i, ] <- sisters[i, order(sisters[i, ])]
	}
	sisters <- unique(sisters)
	final <- matrix(ncol = 4, nrow = dim(sisters)[1])
	colnames(final) <- c("overlap", "range_sp1", "range_sp2", "mrca")
	x <- character()
	for (i in 1:dim(sisters)[1]) {
		x[i] <- paste(sisters[i, 1], sisters[i, 2], sep = "_")
	}
	rownames(final) <- x
	rootAge <- max(node.depth.edgelength(tree))
	for (i in 1:dim(sisters)[1]) {
		final[i, ] <- c(overMatrix[sisters[i, 1], sisters[i, 2]], overMatrix[sisters[i, 1], sisters[i, 1]], overMatrix[sisters[i, 
			2], sisters[i, 2]], rootAge - fastHeight(tree, sisters[i, 1], sisters[i, 2]))
	}
	final[is.na(final)] <- 0
	final <- as.data.frame(final)
	degrOver <- numeric()
	for (i in 1:dim(final)[1]) {
		degrOver[i] <- final$overlap[i] * 100/(min(c(final$range_sp1[i], final$range_sp2[i])))
	}
	final <- cbind(final, degrOver)
	final
}
