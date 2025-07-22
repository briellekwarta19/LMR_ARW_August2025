#different inputs for mod_mort (harvest in our case) work
#you can run the code below with different harvest rates
#harv1 (scalar)
#harv2 (vector, each element applied to each patch)
#harv3 (list, with each matrix in list applied to patch and each element within matrix applied to stage)

library(metapopbio)

n_patches <- 2
n_stages <- 2

group_by <- "patches"
(P <-metapopbio::vec.perm(n_stages = n_stages,
                          n_patches = n_patches,
                          group_by = group_by))

# Northern
f11 <- 0.00
f12 <- 0.2556
s11 <- 0.72
s12 <- 0.77

# Southern
f21 <- 0.00
f22 <- 0.1908 
s21 <- 0.72
s22 <- 0.77

# Northern
(B1x <- matrix(c(f11, f12, s11, s12),
               nrow = n_stages,
               ncol = n_stages,
               byrow = TRUE))

# Southern
(B2x <- matrix(c(f21, f22, s21, s22),
               nrow = n_stages,
               ncol = n_stages,
               byrow = TRUE))

(BB <- metapopbio::blk.diag(list(B1x, B2x)))

# Juveniles
dx1 <- 0.27
dx2 <- 1 - dx1

(Mx1 <- matrix(c(dx2, dx1, dx1, dx2),
               nrow = n_patches,
               ncol = n_patches,
               byrow = TRUE))

# Adults
(Mx2 <- diag(nrow = n_patches))

# Block diagonal matrix
(MM <- metapopbio::blk.diag(list(Mx1, Mx2)))

group_by <- "patches"

lh_order <- "move"

(A <-
    metapopbio::spmm.project.matrix(
      P = P,
      BB = BB,
      MM = MM,
      group_by = group_by,
      lh_order = lh_order
    ))



n <- c(
  60, 19,  # Northern patch adults then juveniles
  29, 20   # Southern patch adults then juveniles
)

comment(n) <- group_by  # vector attr for group_by

n_timesteps <- 100
harv_1 <- 0.2
harv_2 <- c(0.0, 0.2)
harv_3 <- list(
  matrix(c(0.2, 0.0), nrow = 1, byrow = T),
  matrix(c(0.2, 0.2), nrow = 1, byrow = T)
)

projs <-
  metapopbio::spmm.project(
    n = n,
    A = A,
    n_timesteps = n_timesteps,
    n_stages = n_stages,
    n_patches = n_patches,
    mod_mort = harv_3,
    BB = BB,
    MM = MM,
    P = P
  )

stage_names <- c("Juv.", "Adults")

patch_names <- c("North", "South")

metapopbio::spmm.plot(
  projections = projs,
  ylabs = "Abundance",
  xlabs = "Years",
  stage_names = stage_names,
  patch_names = patch_names
)


projsmove <-
  metapopbio::spmm.project(
    n = n,
    A = A,
    n_timesteps = n_timesteps,
    n_stages = n_stages,
    n_patches = n_patches,
    mod_mort = harv_3,
    mod_move = 0.2,
    BB = BB,
    MM = MM,
    P = P
  )


metapopbio::spmm.plot(
  projections = projsmove,
  ylabs = "Abundance",
  xlabs = "Years",
  stage_names = stage_names,
  patch_names = patch_names
)


all.equal(c(projs), c(projsmove))



projsmove

