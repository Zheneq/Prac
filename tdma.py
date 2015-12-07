__author__ = 'Root'


def TDMA(a, b, c, f):
    # b - upper diag: [0;n-2]
    # c - diag: [0;n-1]
    # a - lower diag: [1;n-1]

    a, b, c, f = map(lambda k_list: map(float, k_list), (a, b, c, f))

    alpha = [0]
    beta = [0]
    n = len(f)
    x = [0] * n

    for i in xrange(n-1):
        alpha.append(-b[i]/(a[i]*alpha[i] + c[i]))
        beta.append((f[i] - a[i]*beta[i])/(a[i]*alpha[i] + c[i]))

    x[n-1] = (f[n-1] - a[n-1]*beta[n-1])/(c[n-1] + a[n-1]*alpha[n-1])

    for i in reversed(xrange(n-1)):
        x[i] = alpha[i+1]*x[i+1] + beta[i+1]

    return x