# Load libraries
library(parallel)


# Efficient R programming ---------------------------------------------------------

# Useful link: https://csgillespie.github.io/efficientR/programming.html
help.search(pattern = "optimisation|optimization", fields = c("title", "concept"))




# Cores ---------------------------------------------------------------------------

detectCores() # get number of cores


# * Running "options(mc.cores=XXX)", where XXX is the number of
#   threads your computer can run (2 or 4 on most laptops) will make
#   the model fitting *much* faster
options(mc.cores=4)


# Memory and RAM ------------------------------------------------------------------

gc() # get memory summary 

# A rough rule of thumb is that your RAM should be three times the size of your data set.
benchmarkme::get_ram()
memuse::Sys.meminfo()


object.size("Hello World") # get approximate size of object


# brms optimisation -----------------------------------------------------------------

# https://cran.r-project.org/web/packages/brms/vignettes/brms_threading.html