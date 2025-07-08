library(metapopbio)

n_patches <- 2
n_stages <- 2
group_by <- "patches"

(P <-
    metapopbio::vec.perm(n_stages = n_stages,
                         n_patches = n_patches,
                         group_by = group_by))

# Northern
f11 <- 0.00 
f12 <- 0.26
s11 <- 0.72
s12 <- 0.77

# Southern
f21 <- 0.00
f22 <- 0.19  
s21 <- 0.72
s22 <- 0.77

# Northern
(B1x <-
    matrix(c(f11, f12, s11, s12),
           nrow = n_stages,
           ncol = n_stages,
           byrow = TRUE))
# Southern
(B2x <-
    matrix(c(f21, f22, s21, s22),
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
comment(n) <- "patches"  # vector attr for group_by 

n_timesteps <- 100

head(
  projs <-
    metapopbio::spmm.project(
      n = n,
      A = A,
      n_timesteps = n_timesteps,
      n_stages = n_stages,
      n_patches = n_patches
    )
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

metapopbio::spmm.project.matrix.sens(A)
metapopbio::spmm.eig.lambda(A)
metapopbio::spmm.demo.sens(BB, A, P, MM)
metapopbio::spmm.demo.elas(BB, A, P, MM)
metapopbio::spmm.move.sens(MM, A, P, BB)
metapopbio::spmm.move.elas(MM, A, P, BB)

