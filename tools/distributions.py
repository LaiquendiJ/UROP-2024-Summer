import scipy.stats as stats


def zinb_cdf(k, para):
    """
    Calculate the CDF of the Zero-Inflated Negative Binomial (ZINB).

    Returns:
    cdf : float
        The cumulative probability of the ZINB distribution at k.
    """
    
    mu, r, p_zero = para
    p = 1 - mu / (mu+r)
    nb_cdf = stats.nbinom.cdf(k, r, p)

    return 1-p_zero + p_zero * nb_cdf
