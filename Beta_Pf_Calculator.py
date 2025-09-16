
'''
This module will calculate the reliability index (beta) given the probability of failure (Pf) and vice versa.
'''

from scipy.stats import norm

def reliability_index(Pf):
    """
    Calculate the reliability index (beta) given the probability of failure (Pf).
    
    Args:
        Pf (float): Probability of failure (between 0 and 1)
        
    Returns:
        float: Reliability index beta
    """
    if not (0 < Pf < 1):
        raise ValueError("Probability of failure (Pf) must be between 0 and 1 (exclusive).")
    
    beta = -norm.ppf(Pf)
    return beta

def probability_of_failure(beta):
    """
    Calculate the probability of failure (Pf) given the reliability index (beta).
    
    Args:
        beta (float): Reliability index
        
    Returns:
        float: Probability of failure (between 0 and 1)
    """
    Pf = norm.cdf(-beta)
    return Pf