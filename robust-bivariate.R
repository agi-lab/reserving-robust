### Preamble ===================================================================
# Load required packages 

# To install all packages used below, run the commented section of code:
# install.packages(c("robustbase", "aplpack", "MASS", "ChainLadder", "MVN", "sn", "geometry", "ptinpoly", "rgl", "mrfDepth", "aplpack"))

library(robustbase)
library(MASS)
library(ChainLadder)
library(MVN)
library(sn)
library(geometry)
library(ptinpoly)
library(rgl)
library(mrfDepth)
library(aplpack)
library(xtable)

### find_intersect function ====================================================

# This function will find the coordinates of intersection between a hull defined
# by `hull`, and lines between P0 and points in P1

# hull is the convex hull with which we intersect with
# P0 is the median point
# P1 is the xy coordinates of the points we wish to find intersections for

find_intersect <- function(hull, P0, P1) {
    
    # Construct matrix with points that make up the lines of the bag (each row consists of the two vertices of one line)
    lines <- cbind(hull, hull[c(2:nrow(hull), 1), ])
    
    lines.vectors <- lines[ , 3:4] - lines[ , 1:2]
    
    
    lines.normals <- cbind(-lines.vectors[, 2], lines.vectors[, 1])
    
    # Compute the centroid of each lines and check that they lie on the faces (used to obtain outward facing normals)
    
    centroids <- matrix(nrow = nrow(lines), ncol = 2)
    for(q in 1:nrow(lines)){
        centroids[q, ] <- cbind((lines[q, 1] + lines[q, 3]) / 2, 
                                (lines[q, 2] + lines[q, 4]) / 2)
    }
    
    
    # Create small line segments starting at the centroid and extending in the direction of the normals, recording the endpoints
    
    normal.endpoints <- centroids + (1/100000) * lines.normals # 1/100000 is an arbitrary small number
    
    # Check whether the normals are outward facing by observing if the endpoints are inside or outside the hull, reverse the inward facing normals by multiplying by -1
    
    outward.lines.normals <- matrix(nrow = nrow(lines), ncol = 2)
    for (a in 1:nrow(lines)) {
        
        if(pip2d(Vertices = apply(hull, 2, rev), Queries = t(as.matrix(normal.endpoints[a, ]))) == 1){
            outward.lines.normals[a, ] <- -lines.normals[a, ]
        }
        else {
            outward.lines.normals[a, ] <- lines.normals[a, ]
        }
        
    }
    
    
    # Define P0, P1 and Vi for Cyrus-Beck algorithm (AO min, outlying points)
    P0 <- P0
    P1 <- P1
    Vi <- lines[ , 1:2]
    
    # Compute DS, N and D for use in the Cyrus-Beck algorithm
    
    DS <- matrix(nrow = nrow(P1), ncol = 2) # segment vectors between P0 and P1
    for(m in 1:nrow(DS)){
        DS[m, ] <- P1[m, ] - P0
    }
    
    N <- matrix(nrow = nrow(lines), ncol = 1)
    D <- matrix(nrow = nrow(lines), ncol = nrow(P1))
    potential <- matrix(nrow = nrow(lines), ncol = nrow(P1)) 
    
    # potential will be 1 if the line segment connecting the AO min and outlier is not parallel to a given line,
    # 0 if it is parallel and no intersection with that lines should be considered 
    
    for (t in 1:nrow(lines)){
        N[t] <- -(P0 - Vi[t, ]) %*% outward.lines.normals[t, ]  # negative of dot product of  P0 - Vi and normal vector of each lines
        for (u in 1:nrow(P1)) {
            D[t, u] <- DS[u, ] %*% outward.lines.normals[t, ] # dot product of P1 - P0 and normal vector for each point and each lines
            if (D[t, u] != 0){
                potential[t, u] <- 1}
            else {
                potential[t, u] <- 0}
        }
    }
    
    
    # Evaluate tE (value of t where the line segment enters the hull) and tL (value of t where the line segment exits the hull) for each line segment and lines given by N / D
    
    tE <- matrix(0, nrow = nrow(P1)) # initialize at 0, each row represents the tE value of each outlier
    tL <- matrix(1, nrow = nrow(P1)) # initialize at 1, each row represents the tL value of each outlier
    
    # For each outlier, calculate t for each applicable lines, classify as "leaving" or "entering"
    
    for(z in 1:nrow(P1)){
        for(n in which(potential[ , z] == 1)){  # All faces that each outlier segment may intersect
            t <- N[n] / D[n, z] 
            if(D[n, z] < 0){
                tE[z] <- max(tE[z], t)}   # tE value is max of 0 and all t that are "entering" hull
            else { 
                tL[z] <- min(tL[z], t)}   # tL value is min of 1 and all t that are "leaving" hull
        }
    }
    
    
    # Substitute the tL values into parametric equation of line containing P0
    # and P1 to find the intersection points with the hull (adjusted
    # residual values)
    intersect <- matrix(nrow = nrow(P1), ncol = 2)
    for (c in 1:nrow(P1)) {
        intersect[c, ] <- P0 + tL[c] * (P1[c, ] - P0)
    }
    
    return(intersect)
}

### Import data ==================================================================

# Read dataset from the working directory
ShBaMe12.Robust.Bivariate.CL.Data <- read.csv("data/ShBaMe12 Robust Bivariate CL Data.csv")

# Split dataset into relevant parts that we will be using
Cij.1 <- ShBaMe12.Robust.Bivariate.CL.Data$Cij.1
Cij.2 <- ShBaMe12.Robust.Bivariate.CL.Data$Cij.2
Xij.1 <- ShBaMe12.Robust.Bivariate.CL.Data$Claim.1
Xij.2 <- ShBaMe12.Robust.Bivariate.CL.Data$Claim.2 

### Transform data into triangles ==============================================

# For easy handling we can use the chainladder package and transform the data 
# into actual triangles

i <- rep(1:10,10:1); i <- as.factor(i)
j <- sequence(10:1); j <- as.factor(j)
triangle.names <- c("accyear", "devyear", "incurred_incremental_claims")


## Triangle 1
Frame_Xij.1 <- data.frame(cbind(i, j, Xij.1))
colnames(Frame_Xij.1) = triangle.names
Triangle.Xij.1 <- as.triangle(
    Frame_Xij.1
    , origin = "accyear"
    , dev = "devyear"
    , value = "incurred_incremental_claims"
)

## Triangle 2
Frame_Xij.2 <- data.frame(cbind(i, j, Xij.2))
colnames(Frame_Xij.2) = triangle.names
Triangle.Xij.2 <- as.triangle(
    Frame_Xij.2
    , origin = "accyear"
    , dev = "devyear"
    , value = "incurred_incremental_claims"
)

# Construct cumulative triangles for Triangle 1 & 2
Triangle.Cij.1 <- incr2cum(Triangle.Xij.1)
Triangle.Cij.2 <- incr2cum(Triangle.Xij.2)



### Fit robust GLM Triangle 1 ==================================================

# Fit to Triangle 1
robfit1.1 <- glmrob(
    Xij.1 ~ i + j
    , family = poisson
    , control = glmrobMqle.control(tcc = 1.345, maxit = 500)
    , method = "Mqle"
)

# Compute residuals based on first fit
Residuals.1 <-  (Xij.1 - fitted(robfit1.1)) / sqrt(fitted(robfit1.1))

# Compute 75%-quantile of absolute value of residuals to use as the cutoff 
# value for the second fit
cutrobglm.1 <- quantile(abs(Residuals.1), 0.75)


# Second round fitting using new cutoff value
robfit2.1 <- glmrob(
    Xij.1 ~ i + j
    , family = poisson
    , control = glmrobMqle.control(tcc = cutrobglm.1, maxit = 500)
    , method = "Mqle"
)

# Obtain final residuals based on second robust GLM fit
Final.Residuals.1 <- (ShBaMe12.Robust.Bivariate.CL.Data$Claim.1 - fitted(robfit2.1)) / sqrt(fitted(robfit2.1))

### Fit robust GLM Triangle 2 ==================================================

robfit1.2 <- glmrob(
    Xij.2 ~ i + j
    , family = poisson
    , control = glmrobMqle.control(tcc = 1.345, maxit = 500)
    , method = "Mqle"
)

# Compute residuals based on first fit
Residuals.2 <-  (ShBaMe12.Robust.Bivariate.CL.Data$Claim.2-fitted(robfit1.2))/sqrt(fitted(robfit1.2))

# Compute 75%-quantile of absolute value of residuals to use as the cutoff 
# value for the second fit
cutrobglm.2 <- quantile(abs(Residuals.2), 0.75)

# Second round fitting using new cutoff value
robfit2.2 <- glmrob(
    Xij.2 ~ i + j
    , family = poisson
    , control = glmrobMqle.control(tcc = cutrobglm.2, maxit = 500)
    , method = "Mqle"
)

# Obtain final residuals based on second robust GLM fit
Final.Residuals.2 <- (Xij.2 - fitted(robfit2.2)) / sqrt(fitted(robfit2.2))


## Collate final residuals for both triangles into one matrix
Final.Residuals <- cbind(Final.Residuals.1, Final.Residuals.2)


### Pg. 8 Figure 3: Bagplot ====================================================

pdf("images/bivar_mrfDepth_HD-1_ShBaMe12.pdf")

# Create bagplot of our residuals
bagplot1 <- aplpack::bagplot(
    x = Final.Residuals.1
    , y = Final.Residuals.2
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
    , main = "Bagplot"
)

# Important values from the bagplot - points of outliers, and half-space depths
bp.outliers <- bagplot1$pxy.outlier
bp.non.outliers <- rbind(bagplot1$pxy.bag, bagplot1$pxy.outer)
bp.hdepths <- bagplot1$hdepths

# Four of these observations have the lowest halfspace depth of 1 - plot with
# thick dot
hdepth_1_outliers <- which(rownames(bp.outliers) %in% which(bp.hdepths == 1))
points(
    bp.outliers[hdepth_1_outliers, 1]
    , bp.outliers[hdepth_1_outliers, 2]
    , pch = 19
    , cex = 0.75
    , col = 'red'
)

# Two of the outliers have a halfspace depths of 2 and 3 - marked with red crosses
hdepth_2_3 <- which(rownames(bp.outliers) %in% which(bp.hdepths %in% c(2, 3)))
points(
    bp.outliers[hdepth_2_3, 1]
    , bp.outliers[hdepth_2_3, 2]
    , pch = 4
    , col = "red"
    , cex = 3
)

# One non-outlying observation has a halfspace depth of 1 - marked with a purple cross
non_outlier_hdepth_1 <- which(bp.hdepths == 1)
non_outlier_hdepth_1 <- non_outlier_hdepth_1[!non_outlier_hdepth_1 %in% rownames(bp.outliers)]
points(
    Final.Residuals[non_outlier_hdepth_1, 1]
    , Final.Residuals[non_outlier_hdepth_1, 2]
    , pch = 3
    , col = "violet"
    , cex = 3
)

dev.off()


### Pg 9 Figure 4 Tolerance Ellipse Plots ======================================

### Before adjusting:
MCD <- covMcd(Final.Residuals)
mu <- MCD$center
Sigma <- MCD$cov
MD <- mahalanobis(Final.Residuals, center = mu, cov = Sigma, inverted = FALSE)

pdf("images/bivar_MCD_nonadj_ShBaMe12.pdf", width = 10)

# Plot the ellipse (pg 9, Figure 4(a): Tolerance Ellipses Before Adjusting Outliers)
tolEllipsePlot(
    Final.Residuals
    , m.cov = MCD
    , classic = TRUE
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
)

dev.off()

### Adjust the outliers using the bivariate Winsorization technique
# Set cutoff as defined in paper
MD.cutoff <- qchisq(p = 0.975, df = 2)
MD.c <- qchisq(p = 0.95, df = 2)

# Get adjusted residuals
MCD.adjusted.residuals <- Final.Residuals
for(q in which(MD > MD.cutoff)){
    MCD.adjusted.residuals[q, ] <- min(sqrt(MD.c/MD[q]), 1) * Final.Residuals[q, ] 
}

pdf("images/bivar_MCD_adj_ShBaMe12.pdf", width = 10)

# Plot the ellipse (pg 9, Figure 4(b): Tolerance Ellipses After Adjusting Outliers)
tolEllipsePlot(
    MCD.adjusted.residuals
    , m.cov = covMcd(MCD.adjusted.residuals)
    , classic = TRUE
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
)

dev.off()

### Use adjusted residuals to re-calculate reserves using MCL ~~~~~~~~~~~~~~~~~~
# Backtransform residuals (to incremental claims)

adj.Xij.1 <- sqrt(fitted(robfit2.1)) * MCD.adjusted.residuals[ , 1] + fitted(robfit2.1)

adj.Xij.2 <- sqrt(fitted(robfit2.2)) * MCD.adjusted.residuals[ , 2] + fitted(robfit2.2)

# Set up incremental triangles with backtransformed adjusted residuals
Frame.Xij.1.MCD.adj <- data.frame(cbind(i, j, adj.Xij.1))
colnames(Frame.Xij.1.MCD.adj) = triangle.names
Triangle.Xij.1.MCD.adj <- as.triangle(Frame.Xij.1.MCD.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

Frame.Xij.2.MCD.adj <- data.frame(cbind(i, j, adj.Xij.2))
colnames(Frame.Xij.2.MCD.adj) = triangle.names
Triangle.Xij.2.MCD.adj <- as.triangle(Frame.Xij.2.MCD.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.1.MCD.adj <- incr2cum(Triangle.Xij.1.MCD.adj)
Triangle.Cij.2.MCD.adj <- incr2cum(Triangle.Xij.2.MCD.adj)

# Re-perform multivariate CL and observe results
MCL.MCD <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1.MCD.adj
        , Triangle.Cij.2.MCD.adj
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)


### TODO: Pg. 10 Figure 5, 6: Adjusted Outlyingness Bagplots =========================

# Set a seed for consistency, and use the A.O. function
set.seed(1234)
AO <- adjOutlyingness(Final.Residuals, ndir = 1000)
AO.Final.Residuals <- AO$adjout

# Construct matrix with final residuals and corresponding AO
AO.res.summary <- cbind(Final.Residuals, AO.Final.Residuals)

# Compute the cutoff value
AO.cutoff <- AO$cutoff

# Detect, store and count outliers based on AO cutoff values 1 and 2
AO.outliers.index <- which(AO.Final.Residuals > AO.cutoff)
AO.outliers <- AO.res.summary[AO.outliers.index, 1:2]
count.AO.outliers <- nrow(AO.outliers)

# Store non-outliers
AO.non.outliers <- AO.res.summary[-which(AO.Final.Residuals > AO.cutoff), 1:2] 

# Obtain the data point with the lowest AO (analagous to Tukey median)
Final.Residuals.AO.min <- AO.res.summary[which.min(AO.res.summary[, 3]), 1:2]

# Construct the bag - extract 50% of points with lowest AO and construct
# the convex hull
bag.points.AO <- AO.res.summary[which(AO.res.summary[, 3] <= median(AO.res.summary[, 3])), 1:2]
bag.AO <- convhulln(bag.points.AO, options = "FA")

# Construct the loop - construct the convex hull of all points that are
# not declared as outliers 
loop.AO <- convhulln(AO.non.outliers, options = "FA")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pg.10 Figure 5(a): Adjusted-Outlyingness Bagplot
loophull.AO <- chull(loop.AO$p)
baghull.AO <- chull(bag.AO$p)

pdf("images/bivar_AO-Bagplot_approach3.pdf")

# Plot the bagplot and the required points
plot(bag.points.AO, main = "Adjusted Outlyingness Bagplot", xlab = "Triangle 1", ylab = "Triangle 2", xlim = c(-200, 150), ylim = c(-200, 500))
polygon(loop.AO$p[loophull.AO, ], border = "black", lwd = 2, col = "light blue")
polygon(bag.AO$p[baghull.AO, ],border="black",lwd=2,col = "blue")

points(loop.AO$p,pch=20)
points(AO.outliers, col = "red", pch = "*", cex = 2)
points(Final.Residuals.AO.min, pch = 8, col = "red")

dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Construct the fence by expanding the distance between the bag points
# and AO min by 3x, then construct the convex hull
fence.dist.AO <- 3 * cbind(
    bag.points.AO[, 1] - Final.Residuals.AO.min[1]
    , bag.points.AO[, 2] - Final.Residuals.AO.min[2]
)

fence.points.AO <- cbind(
    fence.dist.AO[, 1] + Final.Residuals.AO.min[1]
    , fence.dist.AO[, 2] + Final.Residuals.AO.min[2]
)

fence.AO <- convhulln(fence.points.AO, options = "FA")
fencehull.AO <- chull(fence.AO$p)

# Detect, store and count the number of outliers using the fence-based
# approach (amount of residuals outside the fence)
AO.outliers.fence <- AO.res.summary[which(pip2d(Vertices = fence.AO$p[rev(fencehull.AO), ], Queries = Final.Residuals) == -1), 1:2]
count.AO.outliers.fence <- nrow(AO.outliers.fence)

# Store non-outliers
AO.non.outliers.fence <- AO.res.summary[-which(pip2d(Vertices = fence.AO$p[rev(fencehull.AO),], Queries = Final.Residuals) == -1), 1:2]

# Construct the loop based on the fence approach
loop.AO.fence <- convhulln(AO.non.outliers.fence, options = "FA")
loophull.AO.fence <- chull(loop.AO.fence$p)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pg. 10 Figure 5(b): Adjusted-Outlyingness Bagplot With Fence

# Plot the bagplot with necessary points

pdf("images/bivar_AO-Bagplot_approach1.pdf")

plot(bag.points.AO, main = "Adjusted Outlyingness Bagplot", xlab = "Triangle 1", ylab = "Triangle 2", xlim = c(-200, 150), ylim = c(-200, 500))
polygon(loop.AO.fence$p[loophull.AO.fence, ], border = "black", lwd = 2, col = "light blue")
polygon(bag.AO$p[baghull.AO, ],border="black",lwd=2,col = "blue")

points(loop.AO.fence$p,pch=20)
points(AO.outliers.fence, col = "red", pch = "*", cex = 2)
points(Final.Residuals.AO.min, pch = 8, col = "red")
polygon(fence.AO$p[fencehull.AO, ], border = "red", lwd = 2)

dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Figure 5(c): Adjusted Outlyingness Bagplot using mrfDepth Cut-off Value

AO.mrfDepth <- adjOutl(Final.Residuals, options = list(ndir = 1000))

AO.Final.Residuals.mrfDepth <- AO.mrfDepth$outlyingnessX

# Construct matrix with final residuals and corresponding AO
AO.res.summary.mrfDepth <- cbind(Final.Residuals, AO.Final.Residuals.mrfDepth)

# Compute the cutoff value
AO.cutoff.mrfDepth <- AO.mrfDepth$cutoff

# Detect, store and count outliers based on AO cutoff values 1 and 2
AO.outliers.mrfDepth <- AO.res.summary.mrfDepth[which(AO.Final.Residuals.mrfDepth > AO.cutoff.mrfDepth), 1:2]

# Store non-outliers
AO.non.outliers.mrfDepth <- AO.res.summary.mrfDepth[-which(AO.Final.Residuals.mrfDepth > AO.cutoff.mrfDepth), 1:2] 

# Obtain the data point with the lowest AO (analagous to Tukey median)
Final.Residuals.AO.min.mrfDepth <- AO.res.summary.mrfDepth[which.min(AO.res.summary.mrfDepth[, 3]), 1:2]

# Construct the bag - extract 50% of points with lowest AO and construct
# the convex hull
bag.points.AO.mrfDepth <- AO.res.summary.mrfDepth[which(AO.res.summary.mrfDepth[, 3] <= median(AO.res.summary.mrfDepth[, 3])), 1:2]
bag.AO.mrfDepth <- convhulln(bag.points.AO.mrfDepth, options = "FA")
baghull.AO.mrfDepth <- chull(bag.AO.mrfDepth$p)

# Construct the loop - construct the convex hull of all points that are
# not declared as outliers 
loop.AO.mrfDepth <- convhulln(AO.non.outliers.mrfDepth, options = "FA")
loophull.AO.mrfDepth <- chull(loop.AO.mrfDepth$p)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pg.10 Figure 5(c): Adjusted-Outlyingness Bagplot using mrfDepth Cut-off Value
### 

pdf("images/bivar_AO-Bagplot_approach2.pdf")

plot(bag.points.AO.mrfDepth, main = "Adjusted Outlyingness Bagplot", xlab = "Triangle 1", ylab = "Triangle 2", xlim = c(-300, 200), ylim = c(-200, 500))
polygon(loop.AO.mrfDepth$p[loophull.AO.mrfDepth, ], border = "black", lwd = 2, col = "light blue")
polygon(bag.AO.mrfDepth$p[baghull.AO.mrfDepth, ], border = "black", lwd = 2, col = "blue")

points(loop.AO.mrfDepth$p, pch = 20)

points(Final.Residuals.AO.min.mrfDepth, pch = 8, col = "red")
points(AO.outliers.mrfDepth, col = "red", pch = '*', cex = 2)

dev.off()

### Pg.11 Figure 7: Bagdistance Illustration ===================================

pdf("images/bivar_mrfDepth_bagdistance.pdf")


bagplot1 <- aplpack::bagplot(
    x = Final.Residuals.1
    , y = Final.Residuals.2
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
    , show.looppoints = FALSE
    , show.bagpoints = FALSE
    , show.whiskers = FALSE
)

# Important coordinates in bagplot: outliers and center of graph
outliers1 <- bagplot1$pxy.outlier
Tukey.Median <- bagplot1$center

# Add orange lines from center to each outlier
segments(x0 = outliers1[, 1], y0 = outliers1[, 2], x1 = Tukey.Median[1], y1 = Tukey.Median[2], col = "orange", lwd = 2)

dev.off()

### TODO: Pg 12 Figure 7 Fence, Loop Adjusted Bagplots ===============================
bagplot1 <- aplpack::bagplot(
    x = Final.Residuals.1
    , y = Final.Residuals.2
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
    , main = "Bagplot"
)

# Important coordinates in bagplot:center of graph
Tukey.Median <- bagplot1$center

# Find Fence # 
baghull1 <- bagplot1$hull.bag
loop1 <- bagplot1$hull.loop

hull.fence <- 3 * cbind(baghull1[, 1] - Tukey.Median[1], baghull1[, 2] - Tukey.Median[2])

hull.fence <- cbind(hull.fence[, 1] + Tukey.Median[1], hull.fence[, 2] + Tukey.Median[2])

# Figure 7(a): Bagplot before Adjusting Outliers with Fence Drawn        

pdf("images/bivar_mrfDepth_bagplot_nonadj_fence.pdf")

plot.bagplot(
    bagplot1
    # , main = "Bagplot"
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
    , xlim = c(-200, 200)
)
polygon(hull.fence, border = "green", lwd=2)
points(bp.outliers, pch = '*', cex = 2, col = 'red')

dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Adjust outliers to fence

fence.outliers.index <- which(pip2d(Vertices = apply(hull.fence, 2, rev), Queries = Final.Residuals) == -1)
fence.outliers <- Final.Residuals[fence.outliers.index, ]

# Find points of intersection between our fence, and lines joining the median
# with each outlier
intersect <- find_intersect(hull.fence, Final.Residuals.AO.min, fence.outliers)

# Adjust the residuals by replacing the outlier values with adjusted values
adj.Final.Residuals.Fence <- Final.Residuals
adj.Final.Residuals.Fence[fence.outliers.index, ] <- intersect[, 1:2]


### Use adjusted residuals to re-calculate reserves using MCL ~~~~~~~~~~~~~~~~~~
# Backtransform residuals (to incremental claims)

adj.Xij.1 <- sqrt(fitted(robfit2.1)) * adj.Final.Residuals.Fence[ , 1] + fitted(robfit2.1)

adj.Xij.2 <- sqrt(fitted(robfit2.2)) * adj.Final.Residuals.Fence[ , 2] + fitted(robfit2.2)

# Set up incremental triangles with backtransformed adjusted residuals
Frame.Xij.1.HD.fence.adj <- data.frame(cbind(i, j, adj.Xij.1))
colnames(Frame.Xij.1.HD.fence.adj) = triangle.names
Triangle.Xij.1.HD.fence.adj <- as.triangle(Frame.Xij.1.HD.fence.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

Frame.Xij.2.HD.fence.adj <- data.frame(cbind(i, j, adj.Xij.2))
colnames(Frame.Xij.2.HD.fence.adj) = triangle.names
Triangle.Xij.2.HD.fence.adj <- as.triangle(Frame.Xij.2.HD.fence.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.1.HD.fence.adj <- incr2cum(Triangle.Xij.1.HD.fence.adj)
Triangle.Cij.2.HD.fence.adj <- incr2cum(Triangle.Xij.2.HD.fence.adj)

# Re-perform multivariate CL and observe results
MCL.HD.fence <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1.HD.fence.adj
        , Triangle.Cij.2.HD.fence.adj
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Figure 7(b): Bagplot After Adjusting Outliers to Fence
# 

pdf("images/bivar_mrfDepth_bagplot_adj-fence.pdf")

bagplot_fence_adjust <- aplpack::bagplot(
    adj.Final.Residuals.Fence
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
    , xlim = c(-200, 200)
    , ylim = c(-200, 400)
)

dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Adjust outliers to loop

loop.outliers.index <- which(pip2d(Vertices = apply(loop1, 2, rev), Queries = Final.Residuals) == -1)
loop.outliers <- Final.Residuals[loop.outliers.index, ]

# Find points of intersection between our loop, and lines joining the median
# with each outlier
intersect <- find_intersect(loop1, Final.Residuals.AO.min, fence.outliers)

# Adjust the residuals by replacing the outlier values with adjusted values
adj.Final.Residuals.loop <- Final.Residuals
adj.Final.Residuals.loop[fence.outliers.index, ] <- intersect[, 1:2]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Figure 9: Bagplot After Adjusting Outliers to Loop

pdf("images/bivar_mrfDepth_bagplot_adj-loop.pdf")

bagplot_loop_adjust <- aplpack::bagplot(
    adj.Final.Residuals.loop
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
    , main = "Bagplot"
    , ylim = c(-40, 60)
)

dev.off()

### Use adjusted residuals to re-calculate reserves using MCL ~~~~~~~~~~~~~~~~~~
# Backtransform residuals (to incremental claims)

adj.Xij.1 <- sqrt(fitted(robfit2.1)) * adj.Final.Residuals.loop[ , 1] + fitted(robfit2.1)

adj.Xij.2 <- sqrt(fitted(robfit2.2)) * adj.Final.Residuals.loop[ , 2] + fitted(robfit2.2)

# Set up incremental triangles with backtransformed adjusted residuals
Frame.Xij.1.HD.loop.adj <- data.frame(cbind(i, j, adj.Xij.1))
colnames(Frame.Xij.1.HD.loop.adj) = triangle.names
Triangle.Xij.1.HD.loop.adj <- as.triangle(Frame.Xij.1.HD.loop.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

Frame.Xij.2.HD.loop.adj <- data.frame(cbind(i, j, adj.Xij.2))
colnames(Frame.Xij.2.HD.loop.adj) = triangle.names
Triangle.Xij.2.HD.loop.adj <- as.triangle(Frame.Xij.2.HD.loop.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.1.HD.loop.adj <- incr2cum(Triangle.Xij.1.HD.loop.adj)
Triangle.Cij.2.HD.loop.adj <- incr2cum(Triangle.Xij.2.HD.loop.adj)

# Re-perform multivariate CL and observe results
MCL.HD.loop <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1.HD.loop.adj
        , Triangle.Cij.2.HD.loop.adj
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)






### Pg. 13 Figure 10: Bagdistance Adjusted Bagplots ==============================

### Calculate the bagdistances

# Find points of the bag
bag.points <- bagplot1$pxy.bag
hull.bag <- bagplot1$hull.bag

# Find out which points are our outliers
outliers1 <- bagplot1$pxy.outlier
outliers.index <- which(paste0(Final.Residuals[, 1], Final.Residuals[, 2]) %in% paste0(outliers1[, 1], outliers1[, 2]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Find points of intersection between our bag, and lines joining the median
# with each outlier
intersect <- find_intersect(hull.bag, Final.Residuals.AO.min, outliers1)

# Compute the numerator term in the bd formula (Euclidean distance between the outliers and Tukey median)
len <- nrow(outliers1)
bd.num <- rep(NA, len)
for (k in 1:len) {
    bd.num[k] <- dist(rbind(outliers1[k, ], Tukey.Median), method = "euclidean")[1]
}


# Compute the denominator term in the bd formula (Euclidean distance between the bag intersection points and Tukey median)
bd.denom <- rep(NA, len)
for (k in 1:len) {
    bd.denom[k] <- dist(rbind(intersect[k, ], Tukey.Median), method = "euclidean")[1]
}

# Compute the bd of outliers by dividing the two terms above
bd.outliers <- bd.num / bd.denom

# Define the factor of the bag we wish to adjust outliers back to
f <- sqrt(qchisq(p = 0.99, df = 2))

# Compute the adjustment factor for outliers for Equation (3.2)
bd.adj.factor <- c()
for (l in 1:length(bd.outliers)) {
    
    if (bd.outliers[l] <= f) {
        bd.adj.factor[l] <- 1
    }
    
    else if (
        (f < bd.outliers[l])
        & (bd.outliers[l] <= 1/2 * (2*f + sqrt(4*f + 1) + 1))
    ) {
        bd.adj.factor[l] <- f / bd.outliers[l]
    }
    
    else if (bd.outliers[l] > 1/2 * (2*f + sqrt(4*f + 1) + 1)) {
        bd.adj.factor[l] <- (f + sqrt(bd.outliers[l])) / bd.outliers[l]
    }
}

# Adjust the residuals using Equation 3.2

adj.Final.Residuals.3.2 <- Final.Residuals
for(l in outliers.index){
    adj.Final.Residuals.3.2[l, ] <- bd.adj.factor[which(outliers.index == l)] * (Final.Residuals[l, ] - Tukey.Median) + Tukey.Median
}


# Figure 9(a) Outliers Adjusted According to Equation (3.2)

pdf("images/bivar_bagdistance_adj1.pdf")

plot.bagplot(
    bagplot1
    , main = "Bagplot"
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
    , xlim = c(-200, 200)
    , show.outlier = FALSE
)
polygon(hull.fence, border="green")

points(
    adj.Final.Residuals.3.2[outliers.index, ]
    , pch = 19
    , col = "orange"
    , cex = 0.75
)

dev.off()

### Use adjusted residuals to re-calculate reserves using MCL ~~~~~~~~~~~~~~~~~~
# Backtransform residuals (to incremental claims)

adj.Xij.1 <- sqrt(fitted(robfit2.1)) * adj.Final.Residuals.3.2[ , 1] + fitted(robfit2.1)

adj.Xij.2 <- sqrt(fitted(robfit2.2)) * adj.Final.Residuals.3.2[ , 2] + fitted(robfit2.2)

# Set up incremental triangles with backtransformed adjusted residuals
i <- rep(1:10,10:1); i <- as.factor(i)
j <- sequence(10:1); j <- as.factor(j)

Frame.Xij.1.bd.3.2.adj <- data.frame(cbind(i, j, adj.Xij.1))
colnames(Frame.Xij.1.bd.3.2.adj) = triangle.names
Triangle.Xij.1.bd.3.2.adj <- as.triangle(Frame.Xij.1.bd.3.2.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

Frame.Xij.2.bd.3.2.adj <- data.frame(cbind(i, j, adj.Xij.2))
colnames(Frame.Xij.2.bd.3.2.adj) = triangle.names
Triangle.Xij.2.bd.3.2.adj <- as.triangle(Frame.Xij.2.bd.3.2.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.1.bd.3.2.adj <- incr2cum(Triangle.Xij.1.bd.3.2.adj)
Triangle.Cij.2.bd.3.2.adj <- incr2cum(Triangle.Xij.2.bd.3.2.adj)

# Re-perform multivariate CL and observe results
MCL.bd.3.2 <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1.bd.3.2.adj
        , Triangle.Cij.2.bd.3.2.adj
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)


# Figure 10 (b) Outliers Adjusted According to Equation (3.3) ~~~~~~~~~~~~~~~~~~

# Set u = 15 as per the paper
u = 15

# Compute the adjustment factor for outliers for Equation (3.3)
bd.adj.factor <- c()
for (l in 1:length(bd.outliers)) {
    
    if (bd.outliers[l] <= f) {
        bd.adj.factor[l] <- 1
    }
    
    else if (
        (
            (f < bd.outliers[l])
            & (bd.outliers[l] <= 1/2 * (2*f + sqrt(4*f + 1) + 1))
        )
        | bd.outliers[l] > u
    ) {
        bd.adj.factor[l] <- f / bd.outliers[l]
    }
    
    else if (bd.outliers[l] > 1/2 * (2*f + sqrt(4*f + 1) + 1)) {
        bd.adj.factor[l] <- (f + sqrt(bd.outliers[l])) / bd.outliers[l]
    }
}

# Adjust the residuals using Equation 3.3

adj.Final.Residuals.3.3 <- Final.Residuals
for(l in outliers.index){
    adj.Final.Residuals.3.3[l, ] <- bd.adj.factor[which(outliers.index == l)] * (Final.Residuals[l, ] - Tukey.Median) + Tukey.Median
}


# Figure 9(b) Outliers Adjusted According to Equation (3.3)

pdf("images/bivar_bagdistance_adj2.pdf")

plot.bagplot(
    bagplot1
    , main = "Bagplot"
    , xlab = "Triangle 1"
    , ylab = "Triangle 2"
    , xlim = c(-200, 200)
    , show.outlier = FALSE
)
polygon(hull.fence, border="green")

points(
    adj.Final.Residuals.3.3[outliers.index, ]
    , pch = 19
    , col = "orange"
    , cex = 0.75
)

dev.off()

### Use adjusted residuals to re-calculate reserves using MCL ~~~~~~~~~~~~~~~~~~
# Backtransform residuals (to incremental claims)

adj.Xij.1 <- sqrt(fitted(robfit2.1)) * adj.Final.Residuals.3.3[ , 1] + fitted(robfit2.1)

adj.Xij.2 <- sqrt(fitted(robfit2.2)) * adj.Final.Residuals.3.3[ , 2] + fitted(robfit2.2)

# Set up incremental triangles with backtransformed adjusted residuals
i <- rep(1:10,10:1); i <- as.factor(i)
j <- sequence(10:1); j <- as.factor(j)

Frame.Xij.1.bd.3.3.adj <- data.frame(cbind(i, j, adj.Xij.1))
colnames(Frame.Xij.1.bd.3.3.adj) = triangle.names
Triangle.Xij.1.bd.3.3.adj <- as.triangle(Frame.Xij.1.bd.3.3.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

Frame.Xij.2.bd.3.3.adj <- data.frame(cbind(i, j, adj.Xij.2))
colnames(Frame.Xij.2.bd.3.3.adj) = triangle.names
Triangle.Xij.2.bd.3.3.adj <- as.triangle(Frame.Xij.2.bd.3.3.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.1.bd.3.3.adj <- incr2cum(Triangle.Xij.1.bd.3.3.adj)
Triangle.Cij.2.bd.3.3.adj <- incr2cum(Triangle.Xij.2.bd.3.3.adj)

# Re-perform multivariate CL and observe results
MCL.bd.3.3 <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1.bd.3.3.adj
        , Triangle.Cij.2.bd.3.3.adj
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)




### TODO: Pg. 11 Table 1: Outlier Detection Results ==================================

# Find halfspace depths of outlier points
bagplot.outliers.index <- which(paste0(Final.Residuals[, 1], Final.Residuals[, 2]) %in% paste0(outliers1[, 1], outliers1[, 2]))
outliers.halfspacedepth <- bagplot1$hdepths[bagplot.outliers.index]

# Find other points with halfspace depth of 1, but weren't classified
# as outliers
other_hd_1.index <- setdiff(which(bagplot1$hdepths == 1), bagplot.outliers.index)
other_hd_1.halfspacedepth <- bagplot1$hdepths[other_hd_1.index]

# Find the MCD Mahalanobis distance for the same points as above
outliers.MCD <- MD[bagplot.outliers.index]
other_hd_1.MCD <- MD[other_hd_1.index]

# Find the AO for the same points as above
outliers.AO <- AO$adjout[bagplot.outliers.index]
other_hd_1.AO <- AO$adjout[other_hd_1.index]

# Find the bagdistance for the same points as above
outliers.bd <- bd.outliers

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



### Calculate the bagdistance for the extra point

# Find points of the bag
bag.points <- bagplot1$pxy.bag
hull.bag <- bagplot1$hull.bag

other_hd_1 <- t(as.matrix(bagplot1$xy[other_hd_1.index, ]))

# Find points of intersection between our bag, and lines joining the median
# with each outlier
intersect <- find_intersect(hull.bag, Final.Residuals.AO.min, outliers1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute the numerator term in the bd formula (Euclidean distance between the outliers and Tukey median)
len <- nrow(other_hd_1)
bd.num <- rep(NA, len)
for (k in 1:len) {
    bd.num[k] <- dist(rbind(other_hd_1[k, ], Tukey.Median), method = "euclidean")[1]
}


# Compute the denominator term in the bd formula (Euclidean distance between the bag intersection points and Tukey median)
bd.denom <- rep(NA, len)
for (k in 1:len) {
    bd.denom[k] <- dist(rbind(intersect[k, ], Tukey.Median), method = "euclidean")[1]
}

# Compute the bd of outliers by dividing the two terms above
other_hd_1.bd <- bd.num / bd.denom


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Finally, compile it all together into a table
table1 <- cbind(
    c(outliers.halfspacedepth, other_hd_1.halfspacedepth)
    , c(outliers.MCD, other_hd_1.MCD)
    , c(outliers.AO, other_hd_1.AO)
    , c(outliers.bd, other_hd_1.bd)
)

# Name the rows
outlier_names <- c(
    "X (6,5)"
    , "X (7,3)"
    , "X (7,4)"
    , "X (8,2)"
    , "X (8,3)"
    , "X (9,1)"
    , "X (9,2)"
)

rownames(table1) <- outlier_names

# Name the columns
colnames(table1) <- c(
    "Bagplot"
    , "MCD"
    , "AO"
    , "bagdistance"
)

# Round the figures in the table and display to console
# 
table1 <- round(table1, digits = 4)
table1

### Extra: AO-based adjustments ================================================

# Find points of intersection between the AO loop, and lines connecting
# the median and outliers
intersect <- find_intersect(loop.AO$hull, Final.Residuals.AO.min, AO.outliers)

# Adjust the residuals by replacing the outlier values with adjusted values above
AO.loop.adjusted.residuals <- Final.Residuals
AO.loop.adjusted.residuals[AO.outliers.index, ] <- intersect[ , 1:2]

### Use adjusted residuals to re-calculate reserves using MCL ~~~~~~~~~~~~~~~~~~
# Backtransform residuals (to incremental claims)

adj.Xij.1 <- sqrt(fitted(robfit2.1)) * AO.loop.adjusted.residuals[ , 1] + fitted(robfit2.1)

adj.Xij.2 <- sqrt(fitted(robfit2.2)) * AO.loop.adjusted.residuals[ , 2] + fitted(robfit2.2)

# Set up incremental triangles with backtransformed adjusted residuals
i <- rep(1:10,10:1); i <- as.factor(i)
j <- sequence(10:1); j <- as.factor(j)
Frame.Xij.1.AO.loop.adj <- data.frame(cbind(i, j, adj.Xij.1))
colnames(Frame.Xij.1.AO.loop.adj) = triangle.names
Triangle.Xij.1.AO.loop.adj <- as.triangle(Frame.Xij.1.AO.loop.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

Frame.Xij.2.AO.loop.adj <- data.frame(cbind(i, j, adj.Xij.2))
colnames(Frame.Xij.2.AO.loop.adj) = triangle.names
Triangle.Xij.2.AO.loop.adj <- as.triangle(Frame.Xij.2.AO.loop.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.1.AO.loop.adj <- incr2cum(Triangle.Xij.1.AO.loop.adj)
Triangle.Cij.2.AO.loop.adj <- incr2cum(Triangle.Xij.2.AO.loop.adj)

# Re-perform multivariate CL and observe results
MCL.AO.loop <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1.AO.loop.adj
        , Triangle.Cij.2.AO.loop.adj
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)






# Find points of intersection between the AO fence, and lines connecting
# the median and outliers
intersect <- find_intersect(fence.AO$hull, Final.Residuals.AO.min, AO.outliers)

# Adjust the residuals by replacing the outlier values with adjusted values above
AO.fence.adjusted.residuals <- Final.Residuals
AO.fence.adjusted.residuals[AO.outliers.index, ] <- intersect[ , 1:2]

### Use adjusted residuals to re-calculate reserves using MCL ~~~~~~~~~~~~~~~~~~
# Backtransform residuals (to incremental claims)

adj.Xij.1 <- sqrt(fitted(robfit2.1)) * AO.fence.adjusted.residuals[ , 1] + fitted(robfit2.1)

adj.Xij.2 <- sqrt(fitted(robfit2.2)) * AO.fence.adjusted.residuals[ , 2] + fitted(robfit2.2)

# Set up incremental triangles with backtransformed adjusted residuals
Frame.Xij.1.AO.fence.adj <- data.frame(cbind(i, j, adj.Xij.1))
colnames(Frame.Xij.1.AO.fence.adj) = triangle.names
Triangle.Xij.1.AO.fence.adj <- as.triangle(Frame.Xij.1.AO.fence.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

Frame.Xij.2.AO.fence.adj <- data.frame(cbind(i, j, adj.Xij.2))
colnames(Frame.Xij.2.AO.fence.adj) = triangle.names
Triangle.Xij.2.AO.fence.adj <- as.triangle(Frame.Xij.2.AO.fence.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.1.AO.fence.adj <- incr2cum(Triangle.Xij.1.AO.fence.adj)
Triangle.Cij.2.AO.fence.adj <- incr2cum(Triangle.Xij.2.AO.fence.adj)

# Re-perform multivariate CL and observe results
MCL.AO.fence <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1.AO.fence.adj
        , Triangle.Cij.2.AO.fence.adj
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)


# Find points of intersection between the AO fence, and lines connecting
# the median and outliers
intersect <- find_intersect(loop.AO.mrfDepth$hull, Final.Residuals.AO.min, AO.outliers)

# Adjust the residuals by replacing the outlier values with adjusted values above
AO.mrfDepth.loop.adjusted.residuals <- Final.Residuals
AO.mrfDepth.loop.adjusted.residuals[AO.outliers.index, ] <- intersect[ , 1:2]

### Use adjusted residuals to re-calculate reserves using MCL ~~~~~~~~~~~~~~~~~~
# Backtransform residuals (to incremental claims)

adj.Xij.1 <- sqrt(fitted(robfit2.1)) * AO.mrfDepth.loop.adjusted.residuals[ , 1] + fitted(robfit2.1)

adj.Xij.2 <- sqrt(fitted(robfit2.2)) * AO.mrfDepth.loop.adjusted.residuals[ , 2] + fitted(robfit2.2)

# Set up incremental triangles with backtransformed adjusted residuals
Frame.Xij.1.AO.mrfDepth.loop.adj <- data.frame(cbind(i, j, adj.Xij.1))
colnames(Frame.Xij.1.AO.mrfDepth.loop.adj) = triangle.names
Triangle.Xij.1.AO.mrfDepth.loop.adj <- as.triangle(Frame.Xij.1.AO.mrfDepth.loop.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

Frame.Xij.2.AO.mrfDepth.loop.adj <- data.frame(cbind(i, j, adj.Xij.2))
colnames(Frame.Xij.2.AO.mrfDepth.loop.adj) = triangle.names
Triangle.Xij.2.AO.mrfDepth.loop.adj <- as.triangle(Frame.Xij.2.AO.mrfDepth.loop.adj, origin = "accyear", dev = "devyear", value = "incurred_incremental_claims")

# Construct cumulative triangles from adjusted incremental triangles
Triangle.Cij.1.AO.mrfDepth.loop.adj <- incr2cum(Triangle.Xij.1.AO.mrfDepth.loop.adj)
Triangle.Cij.2.AO.mrfDepth.loop.adj <- incr2cum(Triangle.Xij.2.AO.mrfDepth.loop.adj)

# Re-perform multivariate CL and observe results
MCL.AO.mrfDepth.loop <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1.AO.mrfDepth.loop.adj
        , Triangle.Cij.2.AO.mrfDepth.loop.adj
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)

### TODO: Pg. 14 Result Tables =======================================================

# Find the indices of the points relevant to our table
table_points_indices <- c(bagplot.outliers.index, other_hd_1.index)
table_colnames <- c(
    "Initial"
    , "MCD"
    , "BP-Fence"
    , "BP-Loop"
    , "AO-Fence"
    , "AO-Loop"
    , "AO-mrfDepth"
    , "bd (no limit)"
    , "bd (limit)"
)

# Get "MCD" column
outliers.MCD.1.adj <- Frame.Xij.1.MCD.adj[table_points_indices, 3]
outliers.MCD.2.adj <- Frame.Xij.2.MCD.adj[table_points_indices, 3]

# Get "Bagplot-Fence" and "Bagplot-Loop" column
outliers.HD.fence.1.adj <- Frame.Xij.1.HD.fence.adj[table_points_indices, 3]
outliers.HD.fence.2.adj <- Frame.Xij.2.HD.fence.adj[table_points_indices, 3]

outliers.HD.loop.1.adj <- Frame.Xij.1.HD.loop.adj[table_points_indices, 3]
outliers.HD.loop.2.adj <- Frame.Xij.2.HD.loop.adj[table_points_indices, 3]

# Get "AO-Fence" and "AO-Loop" column
outliers.AO.loop.1.adj <- Frame.Xij.1.AO.loop.adj[table_points_indices, 3]
outliers.AO.loop.2.adj <- Frame.Xij.2.AO.loop.adj[table_points_indices, 3]

outliers.AO.fence.1.adj <- Frame.Xij.1.AO.fence.adj[table_points_indices, 3]
outliers.AO.fence.2.adj <- Frame.Xij.2.AO.fence.adj[table_points_indices, 3]

# Get "AO-mrfDepth" column
outliers.AO.mrfDepth.loop.1.adj <- Frame.Xij.1.AO.mrfDepth.loop.adj[table_points_indices, 3]
outliers.AO.mrfDepth.loop.2.adj <- Frame.Xij.2.AO.mrfDepth.loop.adj[table_points_indices, 3]

# Get "bd (no limit)" column
outliers.bd.nolimit.1.adj <- Frame.Xij.1.bd.3.2.adj[table_points_indices, 3]
outliers.bd.nolimit.2.adj <- Frame.Xij.2.bd.3.2.adj[table_points_indices, 3]

# Get "bd (limit)" column
outliers.bd.limit.1.adj <- Frame.Xij.1.bd.3.3.adj[table_points_indices, 3]
outliers.bd.limit.2.adj <- Frame.Xij.2.bd.3.3.adj[table_points_indices, 3]

### Table 2: Triangle 1 Outlier Adjustment Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get "Initial" column
initial_triangle1 <- Xij.1[table_points_indices]

# Compile all together into Table 2
table2 <- cbind(
    initial_triangle1
    , outliers.MCD.1.adj
    , outliers.HD.fence.1.adj
    , outliers.HD.loop.1.adj
    , outliers.AO.fence.1.adj
    , outliers.AO.loop.1.adj
    , outliers.AO.mrfDepth.loop.1.adj
    , outliers.bd.nolimit.1.adj
    , outliers.bd.limit.1.adj
)

colnames(table2) <- table_colnames
rownames(table2) <- outlier_names

# For any entry in our table where the value has not changed from the initial,
# set it as NA since it is not an adjustment then
for (i in 1:nrow(table2)) {
    for (j in 2:ncol(table2)) {
        if (table2[i, j] == table2[i, 1]) {
            table2[i, j] <- NA
        }
    }
}

# Round figures in table and output to console
# 
table2 <- round(table2, digits = 4)
table2


### Table 3: Triangle 2 Outlier Adjustment Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get "Initial" column
initial_triangle2 <- Xij.2[table_points_indices]

# Combine all together into Table 3
table3 <- cbind(
    initial_triangle2
    , outliers.MCD.2.adj
    , outliers.HD.fence.2.adj
    , outliers.HD.loop.2.adj
    , outliers.AO.fence.2.adj
    , outliers.AO.loop.2.adj
    , outliers.AO.mrfDepth.loop.2.adj
    , outliers.bd.nolimit.2.adj
    , outliers.bd.limit.2.adj
)

colnames(table3) <- table_colnames
rownames(table3) <- outlier_names

# For any entry in our table where the value has not changed from the initial,
# set it as NA since it is not an adjustment then
for (i in 1:nrow(table3)) {
    for (j in 2:ncol(table3)) {
        if (table3[i, j] == table3[i, 1]) {
            table3[i, j] <- NA
        }
    }
}

# Round figures in table and output to console

table3 <- round(table3, digits = 4)
table3

### Table 4: Bivariate Example Reserves ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get values for "Original" row
MCL.Original <- MultiChainLadder2(
    Triangles = list(
        Triangle.Cij.1
        , Triangle.Cij.2
    )
    , model = "MCL"
    , mse.method = "Independence"
    , control = systemfit::systemfit.control(methodResidCov = "Theil", maxiter = 300)
)

# Initialise Table 4
table4 <- matrix(NA, nrow = 9, ncol = 6)

MCLs <- c(
    "MCL.Original"
    , "MCL.MCD"
    , "MCL.HD.fence"
    , "MCL.HD.loop"
    , "MCL.AO.fence"
    , "MCL.AO.loop"
    , "MCL.AO.mrfDepth.loop"
    , "MCL.bd.3.2"
    , "MCL.bd.3.3"
)

# Get values for each of the rows in Table 4 - these values were fit and
# calculated throughout the previous code, in their relevant sections
for (i in 1:3) {
    for (j in 1:length(MCLs)) {
        table4[j, 2 * i - 1] <- summary(get(MCLs[j]))$IBNR[11, i]
        table4[j, 2 * i] <- summary(get(MCLs[j]))$S.E.Ult[11, i]
    }
    
}

colnames(table4) <- c(
    "Triangle 1 Reserve"
    , "Triangle 1 RMSE"
    , "Triangle 2 Reserve"
    , "Triangle 2 RMSE"
    , "Total Reserve"
    , "Total RMSE"
)

rownames(table4) <- c(
    "Original"
    , "MCD"
    , "Bagplot-Fence"
    , "Bagplot-Loop"
    , "AO-Fence"
    , "AO-Loop"
    , "AO-mrfDepth"
    , "bd (no limit)"
    , "bd (limit)"
)

# Round figures in table and output to console
# 

table4 <- round(table4, digits = 0)
table4

# print(xtable(table1, digits = 4), file="tables/bivar_outlier-detection.txt", booktabs=T)
# print(xtable(table2, digits = 0), file="tables/bivar_tri1-outlier-adj.txt", booktabs=T)
# print(xtable(table3, digits = 0), file="tables/bivar_tri2-outlier-adj.txt", booktabs=T)
# print(xtable(table4, digits = 0), file="tables/bivar_reserve-table.txt", booktabs=T)

write.csv(table1, "tables/bivar_outlier-detection.csv")
write.csv(table2, "tables/bivar_tri1-outlier-adj.csv")
write.csv(table3, "tables/bivar_tri2-outlier-adj.csv")
write.csv(table4, "tables/bivar_reserve-table.csv")