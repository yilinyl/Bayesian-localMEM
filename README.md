# Bayesian local exchangeability design for phase II basket trial

## Abstract

We propose an information borrowing strategy for 
the design and monitoring of phase II basket trials based on the local multisource exchangeability assumption between baskets (disease types). In our proposed local-MEM framework,  information borrowing is only allowed to occur locally, i.e., among baskets with similar response rate and the amount of information borrowing is determined by the level of similarity in response rate, whereas baskets not considered similar are not allowed to share information. We construct a two-stage design for phase II basket trials using the proposed strategy. The proposed method is compared to competing Bayesian methods and Simon's two-stage design in a variety of simulation scenarios. We demonstrate the proposed method is able to maintain the family-wise type I error rate at a reasonable level and has desirable basket-wise power compared to Simon's two-stage design.  In addition, our method is computationally efficient compared to existing Bayesian methods in that the posterior profiles of interest can be derived explicitly without the need for sampling algorithms.

## Code Overview

* `basket.R` implements the core idea of our method including partitioning of baskets into blocks, posterior inference under the proposed local information borrowing strategy
* `prior_sensitivity_analysis.R` gives an example of simulation on fixed design. In this example, we perform 5000 simulations on 4 baskets with maximum sample size all set to 19. The response rates are set to be 0.15 under the null (historical controls) and 0.45 under the alternative (target) for all baskets. Three different priors mentioned in our paper are considered here.
* `sim_2_stg_local_mem.R` performs 5000 simulations on a 2-stage design with 4 baskets. The response rates are set to be 0.15 under the null and 0.45 under the alternative for all baskets. The sample size of each basket is 10 for interim analysis and 16 for final stage analysis, which is consistent with Simon's 2-stage minimax design (please find more details in our paper).
