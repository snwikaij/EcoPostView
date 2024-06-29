# an Ecological Posterior View<br />
This R-package contains functiosn developed for my own work and personal interest. Considering the amount of time refer either to the subsequently use article for publication. Most (if not all) functions use JAGS. This would 
require you to install it first (https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/). The package is focused on visualisation of the posterior and largely on how the (aquatic) ecologist would look at "things". That being said
any question, positive feedback are welcome as well as suggested improvements which benefit the function of the R-package.

Currently it seems a "trend" to be either a dogmatic Bayesian or Frequentist. I completely distantiate myself from thise things as the flaws of one are the benefits of the other and the benefit of the others is the flaw of the other.
Such dogmatic stands are both are in large sense rather unfortionate and largely due to missinterpretations, nihilismn and dogmatism. It is perfectly reasonable to follow the frequentist view, but one has to use the langauge and understanding of it (I tried to explain it a bit here https://snwikaij.shinyapps.io/shiny/). The frequentist cannot make probability statements about the hypothesis, parameters or models ("effect"). 
Choice of words as rejecting, accepting, by chance or the probability of θ falling within x% of the credibility bounds is Bayesian and is often refered to as abduction (or retroduction). Even Fisher was quite clear about it
but his words seem somehow lost (https://github.com/snwikaij/Data/tree/main/Statistical_inference_and_induction). That being said, the Bayesian attaches probabilities to hypothesis, parameters and models via the prior. Then the
Bayesian can state that P(Model|Data) or using Bayes Factor (BF) compare which model most likely generated the data BF10=(P(Model1|Data, Info)/P(Model0|Data, Info))/(P(Model1|Info)/P(Model2|Info). In any case it has rather useful
and pragmatic implication to use Bayesian statistics for ecology. It can then be used to add credibility to the outcome of the data based on prior information (knowledge).

Objective langauge, which I believe in Sander Pierce words would be termed induction. Relies also on prior assumptions, these assumption are however not introduce into the model as the Bayesian does via the prior. 
The researcher has to create ad proper "random" sampling strategy targeting to obtain samples that are (approximattely) independent and identically distributed. When observing the data and the statistic likelihood=P(Data|Model) it 
sees these as an infinity set of "likelihoods" mean study: x̄1={xi, ... xn}/n, x̄2={xi, ... xn}/n ..., Infinity from an infinite set of random variables "true" mean: X̄={X, ..., X}/N where N is infinitely large.
In simplicity 1/2=0.5, 1/3=0.33, 1/4=0.25, then 1/infinity=0. This is important because if state with more samples we get closer to the "true" mean (X̄) we would be Bayesian. Assume we have a large vase with 100 glass marbles
50 read (p) and 50 (q) then taking 90 samples would give use a large amount of "confirmation" about what p/q would be. But if the population infinite we cannot say that a large sample would give use more information on X̄. It
does give us information indeed but only about the data given a model as P(Data|Model) or the data relative to a fixed point P(Data>x|Model) (the latter the p-value). Hence, frequentist statistics has important consideration
when one wants to stay objective one needs to apply infinite many studies to conclude. Ofcourse a 30 more-or-less would from a pragmatic point be sufficien to state something about our emprical error and how this data might
serve as information against some null model (Model0). Furthermore, it can be used a clear signal-to-noise indication  increadibly usefull in randomized trail studies and laboratory studies. Yet the drawbacks that apriori
the the conditions needs to be generated for the results to be interpreted with "confidence".

This Readme is to short to adress all these things and I will not go further in depth either.
