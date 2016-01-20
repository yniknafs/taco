import numpy as np


def _interpolate(a, b, fraction):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a) * fraction


def scoreatpercentile(a, per, limit=(), interpolation_method='fraction'):
    """
    Adapted from scipy

    Calculate the score at the given `per` percentile of the sequence `a`.

    For example, the score at `per=50` is the median. If the desired quantile
    lies between two data points, we interpolate between them, according to
    the value of `interpolation`. If the parameter `limit` is provided, it
    should be a tuple (lower, upper) of two values. Values of `a` outside
    this (closed) interval will be ignored.

    The `interpolation_method` parameter supports three values, namely
    `fraction` (default), `lower` and `higher`. Interpolation is done only,
    if the desired quantile lies between two data points `i` and `j`. For
    `fraction`, the result is an interpolated value between `i` and `j`;
    for `lower`, the result is `i`, for `higher` the result is `j`.

    Parameters
    ----------
    a : ndarray
        Values from which to extract score.
    per : scalar
        Percentile at which to extract score.
    limit : tuple, optional
        Tuple of two scalars, the lower and upper limits within which to
        compute the percentile.
    interpolation : {'fraction', 'lower', 'higher'}, optional
        This optional parameter specifies the interpolation method to use,
        when the desired quantile lies between two data points `i` and `j`:

        - fraction: `i + (j - i)*fraction`, where `fraction` is the
                    fractional part of the index surrounded by `i` and `j`.
        -lower: `i`.
        - higher: `j`.

    Returns
    -------
    score : float
        Score at percentile.

    See Also
    --------
    percentileofscore

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.arange(100)
    >>> stats.scoreatpercentile(a, 50)
    49.5
    """
    # TODO: this should be a simple wrapper around a well-written quantile
    # function.  GNU R provides 9 quantile algorithms (!), with differing
    # behaviour at, for example, discontinuities.
    values = np.sort(a, axis=0)
    if limit:
        values = values[(limit[0] <= values) & (values <= limit[1])]
    idx = per / 100. * (values.shape[0] - 1)
    if (idx % 1 == 0):
        score = values[idx]
    else:
        if interpolation_method == 'fraction':
            score = _interpolate(values[int(idx)], values[int(idx) + 1],
                                 idx % 1)
        elif interpolation_method == 'lower':
            score = values[np.floor(idx)]
        elif interpolation_method == 'higher':
            score = values[np.ceil(idx)]
        else:
            raise ValueError("interpolation_method can only be 'fraction', "
                             "'lower' or 'higher'")
    return score
