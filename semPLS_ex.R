require(seminr)
require(semPLS)

#read in data
mobi = mobi

#create model
mobi_mm <- constructs(
  composite("Image",         multi_items("IMAG", 1:5), weights = mode_B),
  composite("Expectation",   multi_items("CUEX", 1:3), weights = regression_weights),
  composite("Quality",       multi_items("PERQ", 1:7), weights = mode_A),
  composite("Value",         multi_items("PERV", 1:2), weights = correlation_weights),
  reflective("Satisfaction", multi_items("CUSA", 1:3)),
  reflective("Complaints",   single_item("CUSCO")),
  higher_composite("HOC", c("Value", "Satisfaction"), orthogonal, mode_A),
  interaction_term(iv = "Image", moderator = "Expectation", method =  orthogonal, weights = mode_A),
  reflective("Loyalty",      multi_items("CUSL", 1:3))
)

