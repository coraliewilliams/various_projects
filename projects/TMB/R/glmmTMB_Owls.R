library("glmmTMB")
library("bbmle")
library("ggplot2") 

## cosmetic 
theme_set(theme_bw()+ theme(panel.spacing=grid::unit(0,"lines")))

# Load Owls data set ----------

# The data, taken from Zuur et al. (2009) and ultimately from Roulin and Bersier (2007),
# quantify the number of negotiations among owlets (owl chicks) in different nests prior 
# to the arrival of a provisioning parent as a function of food treatment (deprived or satiated),
# the sex of the parent, and arrival time. The total number of calls from the nest is recorded,
# along with the total brood size, which is used as an offset to allow the use of a Poisson
# response. 
# 
# Since the same nests are measured repeatedly, the nest is used as a random effect.
# The model can be expressed as a zero-inflated generalized linear mixed model (ZIGLMM).
# summary(Owls)
summary(Owls)


# Data cleaning 
# (1) reorder nests by mean negotiations per chick, for plotting purposes; 
# (2) add log brood size variable (for offset); 
# (3) rename response variable and abbreviate one of the input variables.
Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  NCalls=SiblingNegotiation,
                  FT=FoodTreatment)


####################
# Zero-inflated models -------------------
####################

# Fit zero-inflated poisson model ------------------------
fit_zipoisson <- glmmTMB(NCalls~(FT+ArrivalTime)*SexParent+
                           offset(log(BroodSize))+(1|Nest),
                         data=Owls,
                         ziformula=~1,
                         family=poisson)

summary(fit_zipoisson)

# Fit zero-inflated negative binomial -------------------
# variance = µ(1 + µ/k
fit_zinbinom <- update(fit_zipoisson,family=nbinom2)

# Fit zero-inflated negative binomial -------------------
# variance = ϕµ
fit_zinbinom1 <- update(fit_zipoisson,family=nbinom1)


# Relax assumption that the total number of calls is strictly proportional 
# to brood size (set as offset term):
fit_zinbinom1_bs <- update(fit_zinbinom1, . ~ (FT+ArrivalTime)*SexParent+ BroodSize+(1|Nest))


# Get goodness of fit in terms of AIC for all models
AICtab(fit_zipoisson,fit_zinbinom,fit_zinbinom1,fit_zinbinom1_bs)


####################
# Hurdle models -------------------
####################


# In contrast to zero-inflated models, hurdle models treat zero-count and nonzero outcomes
# as two completely separate categories, rather than treating the zero-count outcomes as a
# mixture of structural and sampling zeros. 
#
# glmmTMB includes truncated Poisson and negative binomial familes and hence can fit hurdle models.

fit_hnbinom1 <- update(fit_zinbinom1_bs,
                       ziformula=~.,
                       data=Owls,
                       family=truncated_nbinom1)

# Get all AIC values
AICtab(fit_zipoisson,fit_zinbinom, fit_zinbinom1,fit_zinbinom1_bs, fit_hnbinom1)



####################
# Dispersion -------------------
####################

# zero inflated negative binomial with dispersion
fit_zinbinom1_disp <- update(fit_zinbinom1_bs,
                            dispformula = ~1)

# hurdle model negative binomial with dispersion 
fit_hnbinom1_disp <- update(fit_hnbinom1,
                            dispformula = ~1)

# hurdle model negative binomial with dispersion fix and random
fit_hnbinom1_disp_rand <- update(fit_hnbinom1, dispformula=~1+(1|Nest))


# Get all AIC value
AICtab(fit_zipoisson,fit_zinbinom, fit_zinbinom1,fit_zinbinom1_bs, fit_hnbinom1,
       fit_zinbinom1_disp,fit_hnbinom1_disp)
