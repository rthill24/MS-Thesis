import Beta_Pf_Calculator


##Below is a test to get beta from a given Pf
""" Num_Failures = 10
Num_Ships = 1000
P_f = Num_Failures / Num_Ships """

""" beta_test = Beta_Pf_Calculator.reliability_index(P_f)
print ("Beta_test = ", beta_test) """


##Below is a test to get Pf from a given beta
beta = 3 #as a test value

Pf_test = Beta_Pf_Calculator.probability_of_failure(beta)
print ("Pf_test = ", Pf_test)