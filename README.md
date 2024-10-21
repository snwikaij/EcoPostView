# EcoPostView: Ecological Posterior View<br />
This R package contains functions developed for my own work and personal interest. For detailed usage and examples, please refer to the accompanying publication(s). Most functions rely on JAGS, which requires installation (available here: [JAGS Download](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/)). After installing JAGS, you'll need to install the R2jags package in R. The focus of this package is on visualizing posterior distributions, with an emphasis on aquatic ecological data analysis.

Any questions, feedback, or suggestions for improvements to enhance the functionality of the package are welcome.

In recent years, there’s been a noticeable trend towards adopting strict Bayesian or Frequentist perspectives. I distance myself from such dogmatism, as each approach has its strengths and limitations.Such dogmatic stances are, in a larger sense, rather unfortunate and perhaps a more balanced stance is required. It is perfectly reasonable to adopt the frequentist view, but one must use its language and understanding (I attempted to explain it somewhat here: https://snwikaij.shinyapps.io/shiny/). Hence, if exactness, objectivity and error-control the goal of a study Frequentist approaches clearly are your way to go. But this requires you to abstain from probability statments on effects, effect-sizes and hypothesis (model parameters). Choices of words such as rejecting, accepting, by chance, or the probability of θ (theta) falling within x% of the credibility bounds are Bayesian and are often referred to as abduction (or retroduction). Even Fisher was quite clear about this, but his words seem somewhat lost (https://github.com/snwikaij/Data/tree/main/Statistical_inference_and_induction). 

Frequentist approaches offer clear performance measures, such as α (alpha), coverage%=100%(1-α), β (beta), M- and S-type errors. However, Bayesian methods excel in, flexibility and incorporating prior information and compressing it into posterior distributions, though they lack fixed parameters like θ in Frequentist methods (for assesment of performance). Performance in Bayesian models is often assessed via metrics like Mean Squared Error (MSE), but there’s no direct equivalent to α or β.

This package primarily employs Bayesian methods, but it is not dogmatically Bayesian. Future versions may incorporate Frequentist tools as well.
