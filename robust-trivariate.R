################ Set-up ###
# Load required packages 

library(robustbase)
library(aplpack)
library(MASS)
library(ChainLadder)
library(plot3D)
library(geometry)
library(rgl)
library(ptinpoly)
library(depth)
library(MVN)
library(mrfDepth)
library(xtable)

################################################################################
############### Define required functions
# Cross product function for two 3d vectors

cross3d <- function(a, b){
  crossprod <- cbind(a[2] * b[3] - a[3] * b[2],
                     a[3] * b[1] - a[1] * b[3],
                     a[1] * b[2] - a[2] * b[1])
  return(crossprod)
}

# Convert incr claims to triangles
i <- as.factor(rep(1:40, 40:1))
j <- as.factor(sequence(40:1))
triangle.names <- c("accyear", "devyear", "incurred_incremental_claims")

claims_to_tri <- function(incr_claims) {

    Frame.Xij <- data.frame(cbind(i, j, incr_claims))
    colnames(Frame.Xij) = triangle.names
    Triangle.Xij <- as.triangle(Frame.Xij, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")
    
    return(Triangle.Xij)
}

# Create Dataframe of Multivariate CL analysis for IBNR and SE

MCL_toDataframe <- function(MCL) {
    IBNR.LOB1 <- summary(MCL)$IBNR[41, 1]
    S.E.LOB1 <- summary(MCL)$S.E.Ult[41, 1]
    
    IBNR.LOB3 <- summary(MCL)$IBNR[41, 2]
    S.E.LOB3 <- summary(MCL)$S.E.Ult[41, 2]
    
    IBNR.LOB2 <- summary(MCL)$IBNR[41, 3]
    S.E.LOB2 <- summary(MCL)$S.E.Ult[41, 3]
    
    IBNR.Total <- summary(MCL)$IBNR[41, 4]
    S.E.Total <- summary(MCL)$S.E.Ult[41, 4]
    
    return (data.frame(IBNR.LOB1 = IBNR.LOB1,
                       S.E.LOB1 = S.E.LOB1,
                       IBNR.LOB12 = IBNR.LOB3,
                       S.E.LOB12 = S.E.LOB3,
                       IBNR.LOB2 = IBNR.LOB2,
                       S.E.LOB2 = S.E.LOB2,
                       IBNR.Total = IBNR.Total,
                       S.E.Total = S.E.Total))
}

# Poisson GLM fitting Verdonck Debruyne 2011
# Fit robust (poisson) GLM on LOB1 incremental claims data based on Verdonck and Debruyne, 2011
# First round fitting using default cutoff value tcc = 1.345
# # Compute residuals based on first fit
# # Compute 75%-quantile of absolute value of residuals to use as the cutoff value for the second fit 
# Second round fitting using new cutoff value
# # Obtain final residuals based on second robust GLM fit

VD_rob_fit <- function(incr_claims) {
    robfit <- glmrob(incr_claims ~ i + j, family = poisson, 
                         method = "Mqle", control = glmrobMqle.control(tcc = 1.345, maxit = 1000)) # Non-integer warnings
    Residuals <- (incr_claims - fitted(robfit)) / sqrt(fitted(robfit))
    cutoff <- quantile(abs(Residuals), 0.75)
    
    # Second round fitting using new cutoff value
    robfit.final <- glmrob(incr_claims ~ i + j, family = poisson,
                               method = "Mqle", control = glmrobMqle.control(tcc = cutoff, maxit = 1000))
    Final.Residuals <- (incr_claims - fitted(robfit.final)) / sqrt(fitted(robfit.final))
    
    return (list(final.fit = robfit.final, final.res = Final.Residuals))
}

# Adjust vertices back to faces using Cyrus-Beck line clipping algorithm ###

CB_adjustment <- function(adjust_from_centroid, adjust_from_points, adjust_from_index, adjust_to_hull, adjust_to_points,
                          tL_start = 1) {

      # Construct matrix with points that make up the triangular faces of the fence (each row consists of the three vertices of one face)
    faces <- cbind(adjust_to_points[adjust_to_hull$hull[ , 1] , ], adjust_to_points[adjust_to_hull$hull[ , 2] , ], adjust_to_points[adjust_to_hull$hull[ , 3] , ])
    
    # Compute two vectors on each triangular face of the loop and take the cross product to obtain surface normals for each face
    face.vectors <- cbind(faces[ , 4:6] - faces[ , 1:3], faces[ , 7:9] - faces[ , 1:3])
    face.normals <- matrix(nrow = nrow(faces), ncol = 3)
    for(r in 1:nrow(face.normals)){
      face.normals[r, ] <- cross3d(face.vectors[r, 1:3], face.vectors[r, 4:6])  # Using cross product function defined at start
    }
    
    # Compute the centroid of each triangular face and check that they lie on the faces (used to obtain outward facing normals)
    centroids <- matrix(nrow = nrow(faces), ncol = 3)
    for(q in 1:nrow(faces)){
      centroids[q, ] <- cbind((faces[q, 1] + faces[q, 4] + faces[q, 7]) / 3, 
                              (faces[q, 2] + faces[q, 5] + faces[q, 8]) / 3,
                              (faces[q, 3] + faces[q, 6] + faces[q, 9]) / 3)
    }
    
    # Create small line segments starting at the centroid and extending in the direction of the normals, recording the endpoints
    normal.endpoints <- centroids + (1/100000) * face.normals # 1/100000 is an arbitrary small number
    # Check whether the normals are outward facing by observing if the endpoints are inside or outside the hull, reverse the inward facing normals by multiplying by -1
    outward.face.normals <- matrix(nrow = nrow(faces), ncol = 3)
    for(a in 1:nrow(faces)){
      if(pip3d(Vertices = adjust_to_points, Faces = adjust_to_hull$hull, Queries = t(as.matrix(normal.endpoints[a, ]))) == -1){
        outward.face.normals[a, ] <- face.normals[a, ]
      }
      else if(pip3d(Vertices = adjust_to_points, Faces = adjust_to_hull$hull, Queries = t(as.matrix(normal.endpoints[a, ]))) == 1){
        outward.face.normals[a, ] <- -face.normals[a, ]
      } 
    }
    
    # Define P0, P1 and Vi for Cyrus-Beck algorithm (centroid, outliers, and single vertices on each face of the loop respectively)
    P0 <- adjust_from_centroid
    P1 <- adjust_from_points
    Vi <- faces[ , 1:3]
    
    # Compute DS, N and D for use in the Cyrus-Beck algorithm
    DS <- matrix(nrow = nrow(P1), ncol = 3) # segment vectors between P0 and P1
    for(m in 1:nrow(DS)){
      DS[m, ] <- P1[m, ] - P0
    }
    N <- matrix(nrow = nrow(faces), ncol = 1)
    D <- matrix(nrow = nrow(faces), ncol = length(adjust_from_index))
    potential <- matrix(nrow = nrow(faces), ncol = length(adjust_from_index)) 
    
    # potential will be 1 if the line segment connecting the Tukey median and outlier is not parallel to a given face,
    # 0 if it is parallel and no intersection with that face should be considered 
    
    for(t in 1:nrow(faces)){
      N[t] <- -(P0 - Vi[t, ]) %*% outward.face.normals[t, ]  # negative of dot product of  P0 - Vi and normal vector of each face
      for(u in 1:length(adjust_from_index)){
        D[t, u] <- DS[u, ] %*% outward.face.normals[t, ] # dot product of P1 - P0 and normal vector for each point and each face
        if(D[t, u] != 0){
          potential[t, u] <- 1}
        else {
          potential[t, u] <- 0}
      }
    }
    
    # Evaluate tE (value of t where the line segment enters the hull) and tL (value of t where the line segment exits the hull) for each line segment and face given by N / D
    tE <- matrix(0, nrow = length(adjust_from_index)) # initialize at 0, each row represents the tE value of each outlier
    tL <- matrix(tL_start, nrow = length(adjust_from_index)) # initialize at 1, each row represents the tL value of each outlier
    # Note that we initialize at 2 because the value will be greater than 1, as the line segment between the tukey median and inner points will not reach the outer polyhedron
    # For each outlier, calculate t for each applicable face, classify as "leaving" or "entering"
    for(z in 1:length(adjust_from_index)){
      for(n in which(potential[ , z] == 1)){  # All faces that each outlier segment may intersect
        t <- N[n]/D[n, z] 
        if(D[n, z] < 0){
          tE[z] <- max(tE[z], t)}   # tE value is max of 0 and all t that are "entering" hull
        else { 
          tL[z] <- min(tL[z], t)}   # tL value is min of 1 and all t that are "leaving" hull
      }
    }
    
    # Substitute the tL values into parametric equation of line containing P0 and P1 to find the intersection points with the hull (adjusted residual values)
    intersect <- matrix(nrow = length(adjust_from_index), ncol = 3)
    for(c in 1:length(adjust_from_index)){
      intersect[c, ] <- P0 + tL[c] * (P1[c, ] - P0)
    }
    
    # Check that the normal vectors are normal to face vectors - dot product of face vectors and normals should be 0
    check.normality <- matrix(nrow = nrow(faces), ncol = 2)
    for(s in 1:nrow(face.normals)){
      check.normality[s, ] <- cbind(face.normals[s, ] %*% face.vectors[s, 1:3], face.normals[s, ] %*% face.vectors[s, 4:6])
    }
    
    return (list(adj_val = intersect[ , 1:3],
                 norm_check = sum(check.normality),
                 face_centroid_check = pip3d(Vertices = adjust_to_points, Faces = adjust_to_hull$hull, Queries = centroids),
                 adj_check = pip3d(Vertices = adjust_to_points, Faces = adjust_to_hull$hull, Queries = intersect)))
}

################################################################################
# Import Claims Data into R - Change to whatever directory is appropriate 
trivariate.data <- read.csv("data/Trivariate Data.csv")

### Chain Ladder on original data (negatives included) ###
# Obtain incremental claim vectors for each LOB from data (used in paper)
Xij.LOB1 <- trivariate.data$LOB1
Xij.LOB2 <- trivariate.data$LOB2
Xij.LOB3 <- trivariate.data$LOB3

################################################################################
# Check for negative incremental claims in all LOBs - Negative incremental claims present in all LOBs 
length(which(Xij.LOB1 < 0))
length(which(Xij.LOB2 < 0))
length(which(Xij.LOB3 < 0))

# Transform data into triangles (setup) 
# LOB1 Triangle 
# LOB2 Triangle 
# LOB12 Triangle 
Triangle.Xij.LOB1 <- claims_to_tri(Xij.LOB1)
Triangle.Xij.LOB2 <- claims_to_tri(Xij.LOB2)
Triangle.Xij.LOB3 <- claims_to_tri(Xij.LOB3)

# Construct cumulative triangles from incremental for all LOBs
Triangle.Cij.LOB1 <- incr2cum(Triangle.Xij.LOB1)
Triangle.Cij.LOB2 <- incr2cum(Triangle.Xij.LOB2)
Triangle.Cij.LOB3 <- incr2cum(Triangle.Xij.LOB3)

# Classic Mack CL analysis (negative incremental claims included), store CL development factors 
f.classic.LOB1.neg <- MackChainLadder(Triangle = Triangle.Cij.LOB1, est.sigma = "Mack")$f
f.classic.LOB2.neg <- MackChainLadder(Triangle = Triangle.Cij.LOB2, est.sigma = "Mack")$f
f.classic.LOB3.neg <- MackChainLadder(Triangle = Triangle.Cij.LOB3, est.sigma = "Mack")$f

# Preliminary Multivariate CL analysis - store IBNR and S.E for three lines and total
MCL <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1, Triangle.Cij.LOB3, Triangle.Cij.LOB2), model = "MCL", mse.method = "Independence",
                         control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
Prelim.results <- MCL_toDataframe(MCL)

################################################################################
### Chain Ladder on original data with negative incremental claims set to 0 ###
# Set negative incremental claims to 0 
Xij.LOB1[which(Xij.LOB1 < 0)] <- 0
Xij.LOB2[which(Xij.LOB2 < 0)] <- 0
Xij.LOB3[which(Xij.LOB3 < 0)] <- 0

# Transform data into triangles (incremental) 
# LOB1 Triangle (negative claims set to 0) 
# LOB2 Triangle (negative claims set to 0) 
# LOB12 Triangle (negative claims set to 0) 
Triangle.Xij.LOB1 <- claims_to_tri(Xij.LOB1)
Triangle.Xij.LOB2 <- claims_to_tri(Xij.LOB2)
Triangle.Xij.LOB3 <- claims_to_tri(Xij.LOB3)

# Construct cumulative triangles from incremental (negative claims set to 0) 
Triangle.Cij.LOB1 <- incr2cum(Triangle.Xij.LOB1)
Triangle.Cij.LOB2 <- incr2cum(Triangle.Xij.LOB2)
Triangle.Cij.LOB3 <- incr2cum(Triangle.Xij.LOB3)

# Reperform Classic CL analysis after setting negative incremental claims to 0, store development factors
f.classic.LOB1 <- MackChainLadder(Triangle = Triangle.Cij.LOB1, est.sigma = "Mack")$f
f.classic.LOB2 <- MackChainLadder(Triangle = Triangle.Cij.LOB2, est.sigma = "Mack")$f
f.classic.LOB3 <- MackChainLadder(Triangle = Triangle.Cij.LOB3, est.sigma = "Mack")$f

# Multivariate CL analysis after setting negative incremental claims to 0 - store IBNR and S.E for three lines and total
MCL <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1, Triangle.Cij.LOB3, Triangle.Cij.LOB2), model = "MCL", mse.method = "Independence",
                         control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
Original.results <- MCL_toDataframe(MCL)

### Analyse the magnitude of the change in CL development factors after setting negative incremental claims to 0 ###
max(abs((f.classic.LOB1.neg - f.classic.LOB1)/f.classic.LOB1.neg))
max(abs((f.classic.LOB2.neg - f.classic.LOB2)/f.classic.LOB2.neg))
max(abs((f.classic.LOB3.neg - f.classic.LOB3)/f.classic.LOB3.neg)) # max at 0.002307319 = 0.231%

################################################################################
# Fit robust (poisson) GLM on LOB1 incremental claims data based on Verdonck and Debruyne, 2011
# First round fitting using default cutoff value tcc = 1.345
# # Compute residuals based on first fit
# # Compute 75%-quantile of absolute value of residuals to use as the cutoff value for the second fit 
# Second round fitting using new cutoff value
# # Obtain final residuals based on second robust GLM fit

### Obtain residuals for LOB1 ###
### Obtain residuals for LOB2 ###
### Obtain residuals for LOB3 ###
final.fit.LOB1 <- VD_rob_fit(Xij.LOB1)
final.fit.LOB2 <- VD_rob_fit(Xij.LOB2)
final.fit.LOB3 <- VD_rob_fit(Xij.LOB3)

robfit.final.LOB1 <- final.fit.LOB1$final.fit
robfit.final.LOB2 <- final.fit.LOB2$final.fit
robfit.final.LOB3 <- final.fit.LOB3$final.fit

Final.Residuals.LOB1 <- final.fit.LOB1$final.res
Final.Residuals.LOB3 <- final.fit.LOB3$final.res
Final.Residuals.LOB2 <- final.fit.LOB2$final.res

### Collate final residuals for all 3 LOBs into one matrix ###
Final.Residuals.All <- cbind(Final.Residuals.LOB1, Final.Residuals.LOB3, Final.Residuals.LOB2)

################################################################################
################################################################################
################################################################################
### Adjusted Outlyingness outlier adjustment technique on residuals (setup) ###
# Compute the adjusted outlyingness of the residuals using 250 x 3 (750) directions

set.seed(1234)
AO <- adjOutlyingness(Final.Residuals.All, ndir = 750)
AO.Final.Residuals.All <- AO$adjout # Changes with each run

# Construct matrix with final residuals for all 3 LOBs and corresponding AO
AO.res.summary <- cbind(Final.Residuals.All, AO.Final.Residuals.All)

# Compute the two cutoff values to be used to detect outliers - 1st based on medcouple formula, second based on chisquared 99th percentile
AO.cutoff.1 <- AO$cutoff # This is equation (2.3)
AO.cutoff.2 <- sqrt(qchisq(p = 0.99, df = 3)) * median(AO.Final.Residuals.All)

# Detect, store and count outliers based on AO cutoff values 1 and 2, store Non-outliers
AO.outliers.1 <- AO.res.summary[which(AO.Final.Residuals.All > AO.cutoff.1), 1:3]
AO.outliers.2 <- AO.res.summary[which(AO.Final.Residuals.All > AO.cutoff.2), 1:3]
AO.non.outliers.1 <- AO.res.summary[-which(AO.Final.Residuals.All > AO.cutoff.1), 1:3] 
AO.non.outliers.2 <- AO.res.summary[-which(AO.Final.Residuals.All > AO.cutoff.2), 1:3] 
count.AO.outliers.1 <- nrow(AO.outliers.1)
count.AO.outliers.2 <- nrow(AO.outliers.2)

# Obtain the data point with the lowest AO (analagous to Tukey median)
Final.Residuals.AO.min <- AO.res.summary[which.min(AO.res.summary[ , 4]), 1:3]

# Construct the bag - extract 50% of points with lowest AO and construct the convex hull
bag.points.AO <- AO.res.summary[which(AO.res.summary[ , 4] <= median(AO.res.summary[ , 4])), 1:3] 
bag.AO <- convhulln(bag.points.AO, options = "FA")

# Construct the loop - construct the convex hull of all points that are not declared as outliers using both cutoff points (respectively, as above)
loop.AO.1 <- convhulln(AO.non.outliers.1, options = "FA")
loop.AO.2 <- convhulln(AO.non.outliers.2, options = "FA")

# Construct the fence by expanding the distance between the bag points and AO min by 3x, then construct the convex hull
fence.dist.AO <- 3 * cbind(bag.points.AO[ , 1] - Final.Residuals.AO.min[1], bag.points.AO[ , 2] - Final.Residuals.AO.min[2], bag.points.AO[ , 3] - Final.Residuals.AO.min[3])
fence.points.AO <- cbind(fence.dist.AO[ , 1] + Final.Residuals.AO.min[1], fence.dist.AO[ , 2] + Final.Residuals.AO.min[2], fence.dist.AO[ , 3] + Final.Residuals.AO.min[3])
fence.AO <- convhulln(fence.points.AO, options = "FA")

# Detect, store and count the number of outliers using the fence-based approach (amount of residuals outside the fence), store Non-outliers
AO.outliers.3 <- AO.res.summary[which((pip3d(Vertices = fence.points.AO, Faces = fence.AO$hull, Queries = Final.Residuals.All) == -1)), 1:3]
count.AO.outliers.3 <- nrow(AO.outliers.3)
AO.non.outliers.3 <- AO.res.summary[-which((pip3d(Vertices = fence.points.AO, Faces = fence.AO$hull, Queries = Final.Residuals.All) == -1)), 1:3]

# Construct the loop based on the fence approach
loop.AO.3 <- convhulln(AO.non.outliers.3, options = "FA")

################################################################################
### Adjust the outliers back to the loop (traditional cutoff, AO.cutoff.1, eqn 2.3) using Cyrus-Beck line clipping algorithm ###
# Define variable for indices of outliers for later use (all points beyond the loop)
# (adjust_from_centroid, adjust_from_points, adjust_from_index, adjust_to_hull, adjust_to_points)
AO.outliers.index <- which(pip3d(Vertices = AO.non.outliers.1, Faces = loop.AO.1$hull, Queries = Final.Residuals.All) == -1)
AO.outliers.adj <- CB_adjustment(Final.Residuals.AO.min, AO.outliers.1, AO.outliers.index, loop.AO.1, AO.non.outliers.1)
adj.Final.Residuals.All <- Final.Residuals.All
adj.Final.Residuals.All[AO.outliers.index, ] <- AO.outliers.adj$adj_val[ , 1:3]


### Use adjusted residuals to re-calculate reserves using MCL ###
# Backtransform residuals for LOB1, LOB12, LOB2 (to incremental claims)
adj.Xij.LOB1 <- sqrt(fitted(robfit.final.LOB1)) * adj.Final.Residuals.All[ , 1] + fitted(robfit.final.LOB1)
adj.Xij.LOB3 <- sqrt(fitted(robfit.final.LOB3)) * adj.Final.Residuals.All[ , 2] + fitted(robfit.final.LOB3)
adj.Xij.LOB2 <- sqrt(fitted(robfit.final.LOB2)) * adj.Final.Residuals.All[ , 3] + fitted(robfit.final.LOB2)

# Set up incremental triangles with backtransformed adjusted residuals
# LOB1 Triangle 
# LOB2 Triangle 
# LOB3 Triangle 

Triangle.Xij.LOB1.adj <- claims_to_tri(adj.Xij.LOB1)
Triangle.Xij.LOB3.adj <- claims_to_tri(adj.Xij.LOB3)
Triangle.Xij.LOB2.adj <- claims_to_tri(adj.Xij.LOB2)

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.LOB1.adj <- incr2cum(Triangle.Xij.LOB1.adj)
Triangle.Cij.LOB2.adj <- incr2cum(Triangle.Xij.LOB2.adj)
Triangle.Cij.LOB3.adj <- incr2cum(Triangle.Xij.LOB3.adj)

# Reperform multivariate CL - store IBNR and S.E for three lines and total

MCL.AO.loop <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1.adj, Triangle.Cij.LOB2.adj, Triangle.Cij.LOB3.adj), model = "MCL", mse.method = "Independence",
                                 control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
AO.loop.results <- MCL_toDataframe(MCL.AO.loop)

################################################################################
### Adjust the outliers back to the loop (fence approach) using Cyrus-Beck line clipping algorithm ###
# Define variable for indices of outliers for later use (all points beyond the fence)

AO.outliers.index <- which(pip3d(Vertices = fence.points.AO, Faces = fence.AO$hull, Queries = Final.Residuals.All) == -1)
AO.outliers.adj <- CB_adjustment(Final.Residuals.AO.min, AO.outliers.3, AO.outliers.index, loop.AO.3, AO.non.outliers.3)
adj.Final.Residuals.All <- Final.Residuals.All
adj.Final.Residuals.All[AO.outliers.index, ] <- AO.outliers.adj$adj_val[ , 1:3]

### Use adjusted residuals to re-calculate reserves using MCL ###
# Backtransform residuals for LOB1, LOB12, LOB2 (to incremental claims)
adj.Xij.LOB1 <- sqrt(fitted(robfit.final.LOB1)) * adj.Final.Residuals.All[ , 1] + fitted(robfit.final.LOB1)
adj.Xij.LOB3 <- sqrt(fitted(robfit.final.LOB3)) * adj.Final.Residuals.All[ , 2] + fitted(robfit.final.LOB3)
adj.Xij.LOB2 <- sqrt(fitted(robfit.final.LOB2)) * adj.Final.Residuals.All[ , 3] + fitted(robfit.final.LOB2)

# Set up incremental triangles with backtransformed adjusted residuals
# LOB1 Triangle 
# LOB2 Triangle 
# LOB3 Triangle 

Triangle.Xij.LOB1.adj <- claims_to_tri(adj.Xij.LOB1)
Triangle.Xij.LOB3.adj <- claims_to_tri(adj.Xij.LOB3)
Triangle.Xij.LOB2.adj <- claims_to_tri(adj.Xij.LOB2)

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.LOB1.adj <- incr2cum(Triangle.Xij.LOB1.adj)
Triangle.Cij.LOB2.adj <- incr2cum(Triangle.Xij.LOB2.adj)
Triangle.Cij.LOB3.adj <- incr2cum(Triangle.Xij.LOB3.adj)

# Reperform multivariate CL and observe results - store IBNR and S.E for three lines and total
MCL.AO.fence <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1.adj, Triangle.Cij.LOB2.adj, Triangle.Cij.LOB3.adj), model = "MCL", mse.method = "Independence",
                                  control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
AO.fence.results <- MCL_toDataframe(MCL.AO.fence)

################################################################################
### Adjust the outliers back to the loop (mrfDepth, AO.cutoff.2) using Cyrus-Beck line clipping algorithm ###
# Define variable for indices of outliers for later use (all points beyond the loop)
AO.outliers.index <- which(pip3d(Vertices = AO.non.outliers.2, Faces = loop.AO.2$hull, Queries = Final.Residuals.All) == -1)
AO.outliers.adj <- CB_adjustment(Final.Residuals.AO.min, AO.outliers.2, AO.outliers.index, loop.AO.2, AO.non.outliers.2)
adj.Final.Residuals.All <- Final.Residuals.All
adj.Final.Residuals.All[AO.outliers.index, ] <- AO.outliers.adj$adj_val[ , 1:3]

### Use adjusted residuals to re-calculate reserves using MCL ###
# Backtransform residuals for LOB1, LOB12, LOB2 (to incremental claims)
adj.Xij.LOB1 <- sqrt(fitted(robfit.final.LOB1)) * adj.Final.Residuals.All[ , 1] + fitted(robfit.final.LOB1)
adj.Xij.LOB3 <- sqrt(fitted(robfit.final.LOB3)) * adj.Final.Residuals.All[ , 2] + fitted(robfit.final.LOB3)
adj.Xij.LOB2 <- sqrt(fitted(robfit.final.LOB2)) * adj.Final.Residuals.All[ , 3] + fitted(robfit.final.LOB2)

# Set up incremental triangles with backtransformed adjusted residuals
# LOB1 Triangle 
# LOB2 Triangle 
# LOB12 Triangle 

Triangle.Xij.LOB1.adj <- claims_to_tri(adj.Xij.LOB1)
Triangle.Xij.LOB3.adj <- claims_to_tri(adj.Xij.LOB3)
Triangle.Xij.LOB2.adj <- claims_to_tri(adj.Xij.LOB2)

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.LOB1.adj <- incr2cum(Triangle.Xij.LOB1.adj)
Triangle.Cij.LOB2.adj <- incr2cum(Triangle.Xij.LOB2.adj)
Triangle.Cij.LOB3.adj <- incr2cum(Triangle.Xij.LOB3.adj)

# Reperform multivariate CL and observe results - store IBNR and S.E for three lines and total
MCL.AO.mrfDepth <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1.adj, Triangle.Cij.LOB2.adj, Triangle.Cij.LOB3.adj), model = "MCL", mse.method = "Independence",
                                     control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
AO.mrfDepth.results <- MCL_toDataframe(MCL.AO.mrfDepth)

################################################################################
################################################################################
### Halfspace Depth outlier adjustment technique on residuals (setup) ###
# Calculate halfspace depth of the residuals and compute the Tukey median (point of highest HD)

HD.Final.Residuals.All <- c()
for(l in 1:nrow(Final.Residuals.All)){
  HD.Final.Residuals.All[l] <- depth(u = Final.Residuals.All[l, ], x = Final.Residuals.All, approx = FALSE) * nrow(Final.Residuals.All)
}
HD.res.summary <- cbind(Final.Residuals.All, HD.Final.Residuals.All)
Tukey.Median <- Final.Residuals.All[which(HD.Final.Residuals.All == max(HD.Final.Residuals.All)), ] # Check there is only one value
Tukey.Median <- colSums(Tukey.Median)

# Compute the number of points greater than or equal to each unique HD value
unique.HD <- sort(unique(HD.Final.Residuals.All))  # Ordered unique HD values
D.s <- matrix(nrow = length(unique(HD.Final.Residuals.All)))
for(m in 1:length(unique(HD.Final.Residuals.All))){
  D.s[m] <- length(which(HD.Final.Residuals.All >= unique.HD[m]))
}
regions.summary <- cbind(unique.HD, D.s)  # Second column is the number of points with HD >= each unique HD value 

# Extract the HD values to satisfy #Dk <= n/2 < #Dk-1
HD.Dk_1 <- regions.summary[max(which(D.s > nrow(Final.Residuals.All)/2)), 1]  # First HD value of region with no. of points > 410
HD.Dk <- regions.summary[min(which(D.s <= nrow(Final.Residuals.All)/2)), 1]   # First HD value of region with no. of points <= 410

# Store the points in Dk (inner hull) and Dk-1 (outer hull) to construct the bag (note we need >= 4 points to form a convex hull, but this is not a problem here)
Dk_1.pts <- HD.res.summary[which(HD.res.summary[ , 4] >= HD.Dk_1), 1:3]
Dk.pts <- HD.res.summary[which(HD.res.summary[ , 4] >= HD.Dk), 1:3]

# Store points that lie exactly on the inner and outer hulls
outer.points <- Dk_1.pts[which(pip3d(Vertices = Dk_1.pts, Faces = convhulln(Dk_1.pts), Queries = Dk_1.pts) == 0), ]
inner.points <- Dk.pts[which(pip3d(Vertices = Dk.pts, Faces = convhulln(Dk.pts), Queries = Dk.pts) == 0), ]

# Construct the inner and outer hulls
outer.hull <- convhulln(outer.points, options = "FA")
inner.hull <- convhulln(inner.points, options = "FA")

### Adjust points on the outer hull (outside of inner hull) to the inner hull using Cyrus-Beck (to shrink outer polyhedron) ###
# Define variable for indices of points to be adjusted for later use (all points on the outer hull outside of inner hull)
outer.adjust.index <- which(pip3d(Vertices = inner.points, Faces = inner.hull$hull, Queries = outer.points) == -1)
outer.adjust.points <- outer.points[outer.adjust.index, 1:3]
outer.adj <- CB_adjustment(Tukey.Median, outer.adjust.points, outer.adjust.index, inner.hull, inner.points)
HD.intersect.outer <- outer.points
HD.intersect.outer[outer.adjust.index, ] <- outer.adj$adj_val[ , 1:3]

### Adjust points on the inner hull (inside of the outer hull) to the outer hull using Cyrus-Beck (to expand inner polyhedron) ###
# Define variable for indices of points to be adjusted for later use (all points on the inner hull inside of outer hull)
inner.adjust.index <- which(pip3d(Vertices = outer.points, Faces = outer.hull$hull, Queries = inner.points) == 1)
inner.adjust.points <- inner.points[inner.adjust.index, 1:3]
inner.adj <- CB_adjustment(Tukey.Median, inner.adjust.points, inner.adjust.index, outer.hull, outer.points,
                           tL_start = 2)
HD.intersect.inner <- inner.points
HD.intersect.inner[inner.adjust.index, ] <- inner.adj$adj_val[ , 1:3]

lambda <- (nrow(Final.Residuals.All)/2 - nrow(Dk.pts)) / (nrow(Dk_1.pts) - nrow(Dk.pts))

# Compute the bag points by interpolating between outer and inner polygon points (Miller et al. 2003) and construct the bag
# As inner and outer polygons are defined by different vertices, the matrices defining the polygon have different shape
# Below method finds the corresponding point on the face of the outer/inner polygon for each of the vertices on the inner/outer polygon
#   This is defined by altered polygon HD.intersect
# The interpolation is then performed between the altered polygon and original polygon
h1 <- (1 - lambda) * inner.points + lambda * HD.intersect.inner 
h2 <- (1 - lambda) * HD.intersect.outer + lambda * outer.points
h <- rbind(h1, h2)

bag.points.HD <- h
bag.HD <- convhulln(h, options = "FA")

# Construct the fence by expanding the distance between the bag points and Tukey median by the fence factor: chisquared (3 DOF) 99th percentile, then form the convex hull
fence.dist.HD <- sqrt(qchisq(p = 0.99, df = 3)) * cbind(bag.points.HD[ , 1] - Tukey.Median[1], bag.points.HD[ , 2] - Tukey.Median[2], bag.points.HD[ , 3] - Tukey.Median[3])
fence.points.HD <- cbind(fence.dist.HD[ , 1] + Tukey.Median[1], fence.dist.HD[ , 2] + Tukey.Median[2], fence.dist.HD[ , 3] + Tukey.Median[3])
fence.HD <- convhulln(fence.points.HD, options = "FA")

# Construct the loop, formed by the convex hull of all non-outlying points (all points within the fence)
loop.points.HD <- Final.Residuals.All[-c(which(pip3d(Vertices = fence.points.HD, Faces = fence.HD$hull, Queries = Final.Residuals.All) == -1)), ]  # all points within the fence
loop.HD <- convhulln(loop.points.HD, options = "FA")

# Store and count outliers based on HD, store non-outliers
HD.non.outliers <- loop.points.HD
HD.outliers <- Final.Residuals.All[c(which(pip3d(Vertices = fence.points.HD, Faces = fence.HD$hull, Queries = Final.Residuals.All) == -1)), ]
count.HD.outliers <- nrow(HD.outliers)

### Adjust the outliers back to the loop using Cyrus-Beck line clipping algorithm ###
# Define variable for indices of outliers for later use (all points beyond the loop, note that this is the same as all points beyond the fence)
HD.outliers.index <- which(pip3d(Vertices = loop.points.HD, Faces = loop.HD$hull, Queries = Final.Residuals.All) == -1)
HD.outliers.adj <- CB_adjustment(Tukey.Median, HD.outliers, HD.outliers.index, loop.HD, loop.points.HD)
adj.Final.Residuals.All <- Final.Residuals.All
adj.Final.Residuals.All[HD.outliers.index, ] <- HD.outliers.adj$adj_val[ , 1:3]

### Use adjusted residuals to re-calculate reserves using MCL ###
# Backtransform residuals for LOB1, LOB12, LOB2 (to incremental claims)

adj.Xij.LOB1 <- sqrt(fitted(robfit.final.LOB1)) * adj.Final.Residuals.All[ , 1] + fitted(robfit.final.LOB1)
adj.Xij.LOB2 <- sqrt(fitted(robfit.final.LOB2)) * adj.Final.Residuals.All[ , 3] + fitted(robfit.final.LOB2)
adj.Xij.LOB3 <- sqrt(fitted(robfit.final.LOB3)) * adj.Final.Residuals.All[ , 2] + fitted(robfit.final.LOB3)

# Set up incremental triangles with backtransformed adjusted residuals
# LOB1 Triangle 
# LOB2 Triangle 
# LOB12 Triangle 

Triangle.Xij.LOB1.adj <- claims_to_tri(adj.Xij.LOB1)
Triangle.Xij.LOB2.adj <- claims_to_tri(adj.Xij.LOB2)
Triangle.Xij.LOB3.adj <- claims_to_tri(adj.Xij.LOB3)

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.LOB1.adj <- incr2cum(Triangle.Xij.LOB1.adj)
Triangle.Cij.LOB2.adj <- incr2cum(Triangle.Xij.LOB2.adj)
Triangle.Cij.LOB3.adj <- incr2cum(Triangle.Xij.LOB3.adj)

# Reperform multivariate CL and observe results - store IBNR and S.E for three lines and total

MCL.HD.loop <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1.adj, Triangle.Cij.LOB2.adj, Triangle.Cij.LOB3.adj), model = "MCL", mse.method = "Independence",
                                 control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
HD.loop.results <- MCL_toDataframe(MCL.HD.loop)

################################################################################
### Adjust the outliers back to the fence using Cyrus-Beck line clipping algorithm ###
# Define variable for indices of outliers for later use (all points beyond the loop, note that this is the same as all points beyond the fence)
HD.outliers.index <- which(pip3d(Vertices = loop.points.HD, Faces = loop.HD$hull, Queries = Final.Residuals.All) == -1)
HD.outliers.adj <- CB_adjustment(Tukey.Median, HD.outliers, HD.outliers.index, fence.HD, fence.points.HD)
adj.Final.Residuals.All <- Final.Residuals.All
adj.Final.Residuals.All[HD.outliers.index, ] <- HD.outliers.adj$adj_val[ , 1:3]

### Use adjusted residuals to re-calculate reserves using MCL ###
# Backtransform residuals for LOB1, LOB12, LOB2 (to incremental claims)

adj.Xij.LOB1 <- sqrt(fitted(robfit.final.LOB1)) * adj.Final.Residuals.All[ , 1] + fitted(robfit.final.LOB1)
adj.Xij.LOB2 <- sqrt(fitted(robfit.final.LOB2)) * adj.Final.Residuals.All[ , 3] + fitted(robfit.final.LOB2)
adj.Xij.LOB3 <- sqrt(fitted(robfit.final.LOB3)) * adj.Final.Residuals.All[ , 2] + fitted(robfit.final.LOB3)

# Set up incremental triangles with backtransformed adjusted residuals
# LOB1 Triangle 
# LOB2 Triangle 
# LOB12 Triangle 

Triangle.Xij.LOB1.adj <- claims_to_tri(adj.Xij.LOB1)
Triangle.Xij.LOB2.adj <- claims_to_tri(adj.Xij.LOB2)
Triangle.Xij.LOB3.adj <- claims_to_tri(adj.Xij.LOB3)

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.LOB1.adj <- incr2cum(Triangle.Xij.LOB1.adj)
Triangle.Cij.LOB2.adj <- incr2cum(Triangle.Xij.LOB2.adj)
Triangle.Cij.LOB3.adj <- incr2cum(Triangle.Xij.LOB3.adj)

# Reperform multivariate CL and observe results - store IBNR and S.E for three lines and total

MCL.HD.fence <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1.adj, Triangle.Cij.LOB2.adj, Triangle.Cij.LOB3.adj), model = "MCL", mse.method = "Independence",
                                  control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
HD.fence.results <- MCL_toDataframe(MCL.HD.fence)

################################################################################
################################################################################
### MCD Mahalanobis Distance outlier adjustment technique on residuals (setup) ###
# Define MD cutoff value (chisquared 97.5% percentile with 3 DF)
MD.cutoff <- qchisq(p = 0.975, df = 3)

### Robust method
# Obtain the sample location vector (mu) and sample scale matrix (sigma) using the MCD technique (robust)
MCD.Final.Residuals.All <- covMcd(Final.Residuals.All)
mu.r <- MCD.Final.Residuals.All$center
sigma.r <- MCD.Final.Residuals.All$cov

# Calculate the squared Mahalanobis Distance using the MCD estimates (robust), construct matrix with residual values and MD
MD.Final.Residuals.All.r <- mahalanobis(Final.Residuals.All, center = mu.r, cov = sigma.r) # mahalanobis returns MD^2
MD.res.summary.r <- cbind(Final.Residuals.All, MD.Final.Residuals.All.r)

# Detect, count and store outliers using the cutoff value, store non-outliers (robust)
MD.outliers.r <- MD.res.summary.r[which(MD.Final.Residuals.All.r > MD.cutoff), 1:3]
MD.non.outliers.r <- MD.res.summary.r[-which(MD.Final.Residuals.All.r > MD.cutoff), 1:3]
count.MD.outliers.r <- nrow(MD.outliers.r)

### Non Robust method
# Obtain the sample location vector (mu) and sample scale matrix (sigma) using non-robust method
mu.nr <- cbind(mean(Final.Residuals.All[ , 1]), mean(Final.Residuals.All[ , 2]), mean(Final.Residuals.All[ , 3]))
sigma.nr <- cov(Final.Residuals.All)

# Calculate the squared Mahalanobis Distance using the non-robust estimates, construct matrix with residual values and MD
MD.Final.Residuals.All.nr <- mahalanobis(Final.Residuals.All, center = mu.nr, cov = sigma.nr)
MD.res.summary.nr <- cbind(Final.Residuals.All, MD.Final.Residuals.All.nr)

# Detect, count and store outliers using the cutoff value, store non-outliers (non-robust)
MD.outliers.nr <- MD.res.summary.nr[which(MD.Final.Residuals.All.nr > MD.cutoff), 1:3]
MD.non.outliers.nr <- MD.res.summary.nr[which(MD.Final.Residuals.All.nr <= MD.cutoff), 1:3] 
count.MD.outliers.nr <- nrow(MD.outliers.nr)


### Adjust the outliers using trivariate Winsorization technique: min(sqrt(c/MD(x)), 1) * x ###
MD.c <- qchisq(p = 0.95, df = 3) # c value in the formula

# Note that MD.Final.Residuals.All.r is the squared MD, so we take the sqrt inside the formula
adj.Final.Residuals.All <- Final.Residuals.All
for(q in which(MD.Final.Residuals.All.r > MD.cutoff)){
  adj.Final.Residuals.All[q, ] <- min(sqrt(MD.c/MD.Final.Residuals.All.r[q]), 1) * adj.Final.Residuals.All[q, ] 
}


### Use adjusted residuals to re-calculate reserves using MCL ###
# Backtransform residuals for LOB1, LOB12, LOB2 (to incremental claims)
adj.Xij.LOB1 <- sqrt(fitted(robfit.final.LOB1)) * adj.Final.Residuals.All[ , 1] + fitted(robfit.final.LOB1)
adj.Xij.LOB2 <- sqrt(fitted(robfit.final.LOB2)) * adj.Final.Residuals.All[ , 3] + fitted(robfit.final.LOB2)
adj.Xij.LOB3 <- sqrt(fitted(robfit.final.LOB3)) * adj.Final.Residuals.All[ , 2] + fitted(robfit.final.LOB3)

# Set up incremental triangles with backtransformed adjusted residuals
# LOB1 Triangle 
# LOB2 Triangle 
# LOB3 Triangle 

Triangle.Xij.LOB1.adj <- claims_to_tri(adj.Xij.LOB1)
Triangle.Xij.LOB2.adj <- claims_to_tri(adj.Xij.LOB2)
Triangle.Xij.LOB3.adj <- claims_to_tri(adj.Xij.LOB3)

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.LOB1.adj <- incr2cum(Triangle.Xij.LOB1.adj)
Triangle.Cij.LOB2.adj <- incr2cum(Triangle.Xij.LOB2.adj)
Triangle.Cij.LOB3.adj <- incr2cum(Triangle.Xij.LOB3.adj)

# Reperform multivariate CL and observe results - store IBNR and S.E for three lines and total

MCL.MCD <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1.adj, Triangle.Cij.LOB2.adj, Triangle.Cij.LOB3.adj), model = "MCL", mse.method = "Independence",
                             control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
MCD.results <- MCL_toDataframe(MCL.MCD)

################################################################################
################################################################################
### Bagdistance outlier adjustment technique on residuals (setup) ###
### Compute intersection points of the rays through the Tukey median and outliers, and the bag (HD, above) using Cyrus-Beck Algorithm ###
# Define variable for indices of outliers for later use (all points beyond the loop, note that this is the same as all points beyond the fence)
HD.outliers.index <- which(pip3d(Vertices = loop.points.HD, Faces = loop.HD$hull, Queries = Final.Residuals.All) == -1)
HD.outliers.adj <- CB_adjustment(Tukey.Median, HD.outliers, HD.outliers.index, bag.HD, bag.points.HD)

### Calculate bagdistance for the outliers ###
# Compute the numerator term in the bd formula (Euclidean distance between the outliers and Tukey median)
bd.num <- rep(NA, count.HD.outliers)
for(k in 1:count.HD.outliers){
  bd.num[k] <- dist(rbind(Final.Residuals.All[HD.outliers.index[k], ], Tukey.Median), method = "euclidean")[1] 
}

# Compute the denominator term in the bd formula (Euclidean distance between the bag intersection points and Tukey median)
bd.denom <- rep(NA, count.HD.outliers)
for(k in 1:count.HD.outliers){
  bd.denom[k] <- dist(rbind(HD.outliers.adj$adj_val[k, ], Tukey.Median), method = "euclidean")[1]
}

# Compute bd of outliers by dividing the two terms above
bd.outliers <- bd.num/bd.denom

### Adjust the outliers using the first specified adjustment function (Equation 3.2, no "u" limit) ###
# Define the factor of the bag we wish to adjust outliers back to

f <- sqrt(qchisq(p = 0.99, df = 3))  

# Compute the adjustment factor for outliers 

bd.adj.factor <- c()
for(l in 1:length(bd.outliers)){
  if(bd.outliers[l] <= f){
    bd.adj.factor[l] <- 1
  }
  else if((f < bd.outliers[l]) & (bd.outliers[l] <= 1/2 * (2*f + sqrt(4*f + 1) + 1))){
    bd.adj.factor[l] <- f/bd.outliers[l]
  }
  else if(bd.outliers[l] > 1/2 * (2*f + sqrt(4*f + 1) + 1)){
    bd.adj.factor[l] <- (f + sqrt(bd.outliers[l])) / bd.outliers[l]
  }
}

# Adjust the residuals using Equation 3.2
adj.Final.Residuals.All <- Final.Residuals.All
for(l in HD.outliers.index){
  adj.Final.Residuals.All[l, ] <- bd.adj.factor[which(HD.outliers.index == l)] * (Final.Residuals.All[l, ] - Tukey.Median) + Tukey.Median
}


### Use adjusted residuals to re-calculate reserves using MCL ###
# Backtransform residuals for LOB1, LOB12, LOB2 (to incremental claims)
adj.Xij.LOB1 <- sqrt(fitted(robfit.final.LOB1)) * adj.Final.Residuals.All[ , 1] + fitted(robfit.final.LOB1)
adj.Xij.LOB3 <- sqrt(fitted(robfit.final.LOB3)) * adj.Final.Residuals.All[ , 2] + fitted(robfit.final.LOB3)
adj.Xij.LOB2 <- sqrt(fitted(robfit.final.LOB2)) * adj.Final.Residuals.All[ , 3] + fitted(robfit.final.LOB2)

# Set up incremental triangles with backtransformed adjusted residuals
# LOB1 Triangle 
# LOB2 Triangle 
# LOB12 Triangle 

Triangle.Xij.LOB1.adj <- claims_to_tri(adj.Xij.LOB1)
Triangle.Xij.LOB3.adj <- claims_to_tri(adj.Xij.LOB3)
Triangle.Xij.LOB2.adj <- claims_to_tri(adj.Xij.LOB2)

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.LOB1.adj <- incr2cum(Triangle.Xij.LOB1.adj)
Triangle.Cij.LOB2.adj <- incr2cum(Triangle.Xij.LOB2.adj)
Triangle.Cij.LOB3.adj <- incr2cum(Triangle.Xij.LOB3.adj)

# Reperform multivariate CL and observe results - store IBNR and S.E for three lines and total

MCL.bd.no.limit <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1.adj, Triangle.Cij.LOB2.adj, Triangle.Cij.LOB3.adj), model = "MCL", mse.method = "Independence",
                                     control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
bd.no.limit.results <- MCL_toDataframe(MCL.bd.no.limit)


################################################################################
### Adjust the outliers using the first specified adjustment function (Equation 3.3, imposing "u" limit) ###
# Define the limit "u" 
u <- 10 
# Compute the adjustment factor for outliers (including bringing outlier back to fence if bd exceeds limit "u")

bd.adj.factor.2 <- c()
for(l in 1:length(bd.outliers)){
  if(bd.outliers[l] <= f){
    bd.adj.factor.2[l] <- 1
  }
  else if((f < bd.outliers[l]) & (bd.outliers[l] <= 1/2 * (2*f + sqrt(4*f + 1) + 1))){
    bd.adj.factor.2[l] <- f/bd.outliers[l]
  }
  else if((1/2 * (2*f + sqrt(4*f + 1) + 1) < bd.outliers[l]) & (bd.outliers[l] <= u)){
    bd.adj.factor.2[l] <- (f + sqrt(bd.outliers[l])) / bd.outliers[l]
  }
  else if(bd.outliers[l] > u){
    bd.adj.factor.2[l] <- f/bd.outliers[l]
  }
}

# Adjust the residuals using Equation 3.3
adj.Final.Residuals.All <- Final.Residuals.All
for(l in HD.outliers.index){
  adj.Final.Residuals.All[l, ] <- bd.adj.factor.2[which(HD.outliers.index == l)] * (Final.Residuals.All[l, ] - Tukey.Median) + Tukey.Median
}

### Use adjusted residuals to re-calculate reserves using MCL ###
# Backtransform residuals for LOB1, LOB12, LOB2 (to incremental claims)
adj.Xij.LOB1 <- sqrt(fitted(robfit.final.LOB1)) * adj.Final.Residuals.All[ , 1] + fitted(robfit.final.LOB1)
adj.Xij.LOB3 <- sqrt(fitted(robfit.final.LOB3)) * adj.Final.Residuals.All[ , 2] + fitted(robfit.final.LOB3)
adj.Xij.LOB2 <- sqrt(fitted(robfit.final.LOB2)) * adj.Final.Residuals.All[ , 3] + fitted(robfit.final.LOB2)

# Set up incremental triangles with backtransformed adjusted residuals
# LOB1 Triangle 
# LOB2 Triangle 
# LOB12 Triangle 

Triangle.Xij.LOB1.adj <- claims_to_tri(adj.Xij.LOB1)
Triangle.Xij.LOB3.adj <- claims_to_tri(adj.Xij.LOB3)
Triangle.Xij.LOB2.adj <- claims_to_tri(adj.Xij.LOB2)

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.LOB1.adj <- incr2cum(Triangle.Xij.LOB1.adj)
Triangle.Cij.LOB2.adj <- incr2cum(Triangle.Xij.LOB2.adj)
Triangle.Cij.LOB3.adj <- incr2cum(Triangle.Xij.LOB3.adj)
# Reperform multivariate CL and observe results - store IBNR and S.E for three lines and total

MCL.bd.limit <- MultiChainLadder2(Triangles = list(Triangle.Cij.LOB1.adj, Triangle.Cij.LOB2.adj, Triangle.Cij.LOB3.adj), model = "MCL", mse.method = "Independence",
                                  control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300))
bd.limit.results <- MCL_toDataframe(MCL.bd.limit)

# Define row and column names for table

table.col.names <- c("Reserve (LOB1 Insurer 1)", "RMSE (LOB1 Insurer 1)", "Reserve (LOB3 Insurer 2)", "RMSE (LOB3 Insurer 2)", 
                     "Reserve (LOB2 Insurer 1)", "RMSE (LOB2 Insurer 1)", "Reserve (Total)", "RMSE (Total)")

table.row.names <- c("Original", "AO Loop", "AO Fence", "AO-mrfDepth", "HD Loop", "HD Fence", "MCD", "bd (no limit)", "bd (limit)")

# Construct table

Reserve.summary.table <- rbind(Original.results, AO.loop.results, AO.fence.results, AO.mrfDepth.results, HD.loop.results, HD.fence.results, MCD.results, bd.no.limit.results, bd.limit.results)
colnames(Reserve.summary.table) <- table.col.names
rownames(Reserve.summary.table) <- table.row.names


View(Reserve.summary.table)
write.csv(Reserve.summary.table, "tables/trivar_reserve-table.csv")

# ################################################################################
# ###### Graphing
# ################################################################################

graph = T

if (graph) {

# ### AO plots ### 

# AO Plot Cutoff 1 bag
open3d()
rgl::rgl.viewpoint(60)
rgl::rgl.light(120, 60)

# plot non-outliers in blue, outliers in red (AO cutoff 1, used for example in the paper)
points3d(AO.non.outliers.1, col = "blue", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)
points3d(AO.outliers.1, col = "red", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)

# identify and store which points form the convex hulls (bag, loop, fence) to plot - need to transpose so the order of points in triangles is correct
bag.AO.surf <- t(bag.AO$hull)
loop.AO.surf <- t(loop.AO.1$hull)
fence.AO.surf <- t(fence.AO$hull)

# plot the bag using indices from convex hull commands 
triangles3d(x = bag.points.AO[bag.AO.surf, 1], y = bag.points.AO[bag.AO.surf, 2] , z = bag.points.AO[bag.AO.surf, 3], col = "yellow", alpha = 0.6)

# show axes and grid
grid3d(c("x", "y", "z"))
decorate3d(xlab = "LOB 1", ylab = "LOB 2", zlab = "LOB 3", box = F, axes = F)
rgl::rgl.snapshot("images/trivar_AO-bagplot.png")

# AO plot cutoff 1 
open3d()
rgl::rgl.viewpoint(60)
rgl::rgl.light(120, 60)

points3d(AO.non.outliers.1, col = "blue", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)
points3d(AO.outliers.1, col = "red", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)

triangles3d(x = bag.points.AO[bag.AO.surf, 1], y = bag.points.AO[bag.AO.surf, 2] , z = bag.points.AO[bag.AO.surf, 3], col = "yellow", alpha = 0.6)
triangles3d(x = AO.non.outliers.1[loop.AO.surf, 1], y = AO.non.outliers.1[loop.AO.surf, 2], z = AO.non.outliers.1[loop.AO.surf, 3], col = "#a53900", alpha = 0.5)

grid3d(c("x", "y", "z"))
decorate3d(xlab = "LOB 1", ylab = "LOB 2", zlab = "LOB 3", box = F, axes = F)
rgl::rgl.snapshot("images/trivar_AO-bagplot-loop_approach3.png")


# AO Plot Cutoff 1
open3d()
rgl::rgl.viewpoint(60)
rgl::rgl.light(120, 60)

points3d(AO.non.outliers.1, col = "blue", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)
points3d(AO.outliers.1, col = "red", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)

triangles3d(x = bag.points.AO[bag.AO.surf, 1], y = bag.points.AO[bag.AO.surf, 2] , z = bag.points.AO[bag.AO.surf, 3], col = "yellow", alpha = 0.6)
triangles3d(x = fence.points.AO[fence.AO.surf, 1], y = fence.points.AO[fence.AO.surf, 2], z = fence.points.AO[fence.AO.surf, 3], col = "light blue", alpha = 0.6) # This is the fence

grid3d(c("x", "y", "z"))
decorate3d(xlab = "LOB 1", ylab = "LOB 2", zlab = "LOB 3", box = F, axes = F)
rgl::rgl.snapshot("images/trivar_AO-bagplot-loop.png")


# AO Plot Cutoff 1
rgl::rgl.snapshot("images/trivar_AO-bagplot_1-3comparison.png")
open3d()
rgl::rgl.viewpoint(60)
rgl::rgl.light(120, 60)

points3d(AO.non.outliers.1, col = "blue", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)
points3d(AO.outliers.1, col = "red", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)

triangles3d(x = bag.points.AO[bag.AO.surf, 1], y = bag.points.AO[bag.AO.surf, 2] , z = bag.points.AO[bag.AO.surf, 3], col = "yellow", alpha = 0.6)
triangles3d(x = fence.points.AO[fence.AO.surf, 1], y = fence.points.AO[fence.AO.surf, 2], z = fence.points.AO[fence.AO.surf, 3], col = "light blue", alpha = 0.6)
triangles3d(x = AO.non.outliers.1[loop.AO.surf, 1], y = AO.non.outliers.1[loop.AO.surf, 2], z = AO.non.outliers.1[loop.AO.surf, 3], col = "#a53900", alpha = 0.6)

grid3d(c("x", "y", "z"))
decorate3d(xlab = "LOB 1", ylab = "LOB 2", zlab = "LOB 3", box = F, axes = F)

### HD plots ###
# set up 3d plot device
open3d()
rgl::rgl.viewpoint(60)
rgl::rgl.light(120, 60)

# plot non-outliers in blue, outliers in red
points3d(HD.non.outliers, col = "blue", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)
points3d(HD.outliers, col = "red", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)

# identify and store which points form the convex hulls (bag, loop, fence) to plot - need to transpose so the order of points in triangles is correct
bag.HD.surf <- t(bag.HD$hull)
loop.HD.surf <- t(loop.HD$hull)
fence.HD.surf <- t(fence.HD$hull)

# plot the bag using indices from convex hull commands 
triangles3d(x = bag.points.HD[bag.HD.surf, 1], y = bag.points.HD[bag.HD.surf, 2] , z = bag.points.HD[bag.HD.surf, 3], col = "yellow", alpha = 0.6)

# show axes and grid
grid3d(c("x", "y", "z"))
decorate3d(xlab = "LOB 1", ylab = "LOB 2", zlab = "LOB 3", box = F, axes = F)
rgl::rgl.snapshot("images/trivar_HD-bag.png")

# plot points with bag, loop and fence 
open3d()
rgl::rgl.viewpoint(60)
rgl::rgl.light(120, 60)

points3d(HD.non.outliers, col = "blue", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)
points3d(HD.outliers, col = "red", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)

triangles3d(x = bag.points.HD[bag.HD.surf, 1], y = bag.points.HD[bag.HD.surf, 2] , z = bag.points.HD[bag.HD.surf, 3], col = "yellow", alpha = 0.6)
triangles3d(x = fence.points.HD[fence.HD.surf, 1], y = fence.points.HD[fence.HD.surf, 2], z = fence.points.HD[fence.HD.surf, 3], col = "light blue", alpha = 0.6)
triangles3d(x = loop.points.HD[loop.HD.surf, 1], y = loop.points.HD[loop.HD.surf, 2], z = loop.points.HD[loop.HD.surf, 3], col = "#a53900", alpha = 0.6)

grid3d(c("x", "y", "z"))
decorate3d(xlab = "LOB 1", ylab = "LOB 2", zlab = "LOB 3", box = F, axes = F)
rgl::rgl.snapshot("images/trivar_HD-bag-loop-fence.png")

### MD plots ###
# MCD MD points with robust tolerance ellipse 
# set up 3d plot device
open3d()
rgl::rgl.viewpoint(60)
rgl::rgl.light(120, 60)

# plot non-outliers in blue, outliers in red (using MD cutoff)
points3d(MD.non.outliers.r, col = "blue", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)
points3d(MD.outliers.r, col = "red", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)

# plot the robust tolerance ellipse (97.5% confidence level)
plot3d(ellipse3d(sigma.r, centre = mu.r, t = sqrt(MD.cutoff)), col = "#a53900", alpha = 0.5, add = T)

# show axes and grid
grid3d(c("x", "y", "z"))
decorate3d(xlab = "LOB 1", ylab = "LOB 2", zlab = "LOB 3", box = F, axes = F)
rgl::rgl.snapshot("images/trivar_MCD-ellipse.png")

# MCD ML points with non-robust and robust tolerance ellipses 
# set up 3d plot device
open3d()
rgl::rgl.viewpoint(60)
rgl::rgl.light(120, 60)

# plot non-outliers in blue, outliers in red (using MD cutoff)
points3d(MD.non.outliers.nr, col = "blue", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)
points3d(MD.outliers.nr, col = "red", alpha = 0.7, aspect = c(1, 1, 0.5), size = 5)

# plot the robust and non-robust tolerance ellipses (97.5% confidence level)
plot3d(ellipse3d(sigma.r, centre = mu.r, t = sqrt(MD.cutoff)), col = "#a53900", alpha = 0.5, add = T)
plot3d(ellipse3d(sigma.nr, centre = mu.nr, t = sqrt(MD.cutoff)), col = "light green", alpha = 0.5, add = T)

# show axes and grid
grid3d(c("x", "y", "z"))
decorate3d(xlab = "LOB 1", ylab = "LOB 2", zlab = "LOB 3", box = F, axes = F)
rgl::rgl.snapshot("images/trivar_MCD-ellipse_robust.png")


}

