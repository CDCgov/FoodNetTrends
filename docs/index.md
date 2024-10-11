# Background

The Foodborne Diseases Active Surveillance Network (FoodNet) monitors illnesses caused by enteric and foodborne pathogens across 10 U.S. sites. FoodNet data is used to track trends in these illnesses and to monitor progress toward federal disease reduction goals.

The original model for analyzing FoodNet data faced limitations, such as sensitivity to single-year aberrations and biases toward more populous sites. To address these issues, this enhanced model (FoodNetTrends) was developed using a Bayesian framework, incorporating thin-plate splines and site-specific interactions.

Key improvements include:
- Treating the year as a continuous variable.
- Including site-specific trends.
- Improved ability to handle uncertainty and noisy data.

[User Guide](user_guide.md)

[FAQ](faq.md)

[Github](https://github.com/CDCgov/FoodNetTrends)
